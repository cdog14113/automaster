/* Auto Master DSP. Plain JS, no dependencies. Works in the browser and in Node (for tests).
   Steps: measure (BS.1770 / EBU R128 loudness, true peak) -> gain to target -> lookahead limiter -> verify. */
(function (root) {
  'use strict';

  // ---- K-weighting (ITU-R BS.1770-4), coefficients derived per sample rate ----
  function kWeightCoefs(fs) {
    // Stage 1: high shelf (head effect)
    var f0 = 1681.974450955533, G = 3.999843853973347, Q = 0.7071752369554196;
    var K = Math.tan(Math.PI * f0 / fs), Vh = Math.pow(10, G / 20), Vb = Math.pow(Vh, 0.4996667741545416);
    var a0 = 1 + K / Q + K * K;
    var s1 = { b0: (Vh + Vb * K / Q + K * K) / a0, b1: 2 * (K * K - Vh) / a0, b2: (Vh - Vb * K / Q + K * K) / a0,
               a1: 2 * (K * K - 1) / a0, a2: (1 - K / Q + K * K) / a0 };
    // Stage 2: RLB high-pass
    f0 = 38.13547087602444; Q = 0.5003270373238773;
    K = Math.tan(Math.PI * f0 / fs);
    a0 = 1 + K / Q + K * K;
    var s2 = { b0: 1, b1: -2, b2: 1, a1: 2 * (K * K - 1) / a0, a2: (1 - K / Q + K * K) / a0 };
    return [s1, s2];
  }

  function biquad(x, c) {
    var y = new Float32Array(x.length), x1 = 0, x2 = 0, y1 = 0, y2 = 0;
    for (var i = 0; i < x.length; i++) {
      var xi = x[i];
      var yi = c.b0 * xi + c.b1 * x1 + c.b2 * x2 - c.a1 * y1 - c.a2 * y2;
      y[i] = yi; x2 = x1; x1 = xi; y2 = y1; y1 = yi;
    }
    return y;
  }

  // Integrated loudness in LUFS with absolute (-70) and relative (-10 LU) gating, 400 ms blocks, 75% overlap.
  function integratedLoudness(chs, fs) {
    var co = kWeightCoefs(fs), n = chs[0].length, filt = [];
    for (var c = 0; c < chs.length; c++) filt.push(biquad(biquad(chs[c], co[0]), co[1]));
    var blk = Math.round(0.4 * fs), step = Math.round(0.1 * fs);
    if (n < blk) { // too short for one block: use the whole signal
      blk = n; step = n;
    }
    var nb = Math.floor((n - blk) / step) + 1, z = [];
    // prefix sums of squares per channel (Float64)
    var pre = filt.map(function (f) { var p = new Float64Array(n + 1); for (var i = 0; i < n; i++) p[i + 1] = p[i] + f[i] * f[i]; return p; });
    for (var b = 0; b < nb; b++) {
      var s = b * step, e = s + blk, sum = 0;
      for (var c2 = 0; c2 < pre.length; c2++) sum += (pre[c2][e] - pre[c2][s]) / blk;
      z.push(sum);
    }
    var toL = function (m) { return -0.691 + 10 * Math.log10(m); };
    var abs = z.filter(function (m) { return m > 0 && toL(m) > -70; });
    if (!abs.length) return -Infinity;
    var mean = abs.reduce(function (a, b) { return a + b; }, 0) / abs.length;
    var relT = toL(mean) - 10;
    var rel = abs.filter(function (m) { return toL(m) > relT; });
    if (!rel.length) return -Infinity;
    return toL(rel.reduce(function (a, b) { return a + b; }, 0) / rel.length);
  }

  // ---- True peak: 4x oversampling with windowed-sinc interpolation ----
  var PH = (function () {
    var half = 6, phases = [];
    for (var p = 1; p < 4; p++) {
      var f = p / 4, taps = [], sum = 0;
      for (var k = -half + 1; k <= half; k++) {
        var t = k - f, sinc = Math.sin(Math.PI * t) / (Math.PI * t);
        var w = 0.5 * (1 + Math.cos(Math.PI * t / half));
        taps.push(sinc * w); sum += sinc * w;
      }
      phases.push(taps.map(function (v) { return v / sum; }));
    }
    return { half: half, phases: phases };
  })();

  // Per-sample peak envelope: max(|x[n]|, |interpolated values between n and n+1|) across channels.
  function peakEnvelope(chs) {
    var n = chs[0].length, env = new Float32Array(n), half = PH.half;
    for (var c = 0; c < chs.length; c++) {
      var x = chs[c];
      for (var i = 0; i < n; i++) {
        var m = Math.abs(x[i]);
        for (var p = 0; p < 3; p++) {
          var taps = PH.phases[p], v = 0, j0 = i - half + 1;
          for (var k = 0; k < taps.length; k++) { var j = j0 + k; if (j >= 0 && j < n) v += taps[k] * x[j]; }
          v = Math.abs(v); if (v > m) m = v;
        }
        if (m > env[i]) env[i] = m;
      }
    }
    return env;
  }
  function truePeakDb(chs) {
    var env = peakEnvelope(chs), m = 0;
    for (var i = 0; i < env.length; i++) if (env[i] > m) m = env[i];
    return m > 0 ? 20 * Math.log10(m) : -Infinity;
  }
  function samplePeakDb(chs) {
    var m = 0; chs.forEach(function (x) { for (var i = 0; i < x.length; i++) { var a = Math.abs(x[i]); if (a > m) m = a; } });
    return m > 0 ? 20 * Math.log10(m) : -Infinity;
  }
  // Samples at (or beyond) full scale, a sign of clipping. Counts runs of 3+ in a row as clipping too.
  function clippedSamples(chs) {
    var cnt = 0;
    chs.forEach(function (x) { for (var i = 0; i < x.length; i++) if (Math.abs(x[i]) >= 0.9999) cnt++; });
    return cnt;
  }

  // ---- Lookahead brick-wall limiter (offline) ----
  // g[n] = gain needed so the peak at n stays under the ceiling. min-filter over the lookahead
  // window, then a moving average of the same length: the result is <= g[n] everywhere (so no peak gets through),
  // with smooth ramps (no distortion from sudden gain steps). Release is slower and smooth.
  function limiterGain(env, ceilLin, fs, lookMs, relMs) {
    var n = env.length, L = Math.max(2, Math.round(lookMs / 1000 * fs));
    var g = new Float32Array(n);
    for (var i = 0; i < n; i++) g[i] = env[i] > ceilLin ? ceilLin / env[i] : 1;
    // sliding min over g[i..i+L-1] (monotonic deque)
    var m = new Float32Array(n), dq = new Int32Array(n), h = 0, t = 0;
    for (var k = n - 1; k >= 0; k--) {
      while (t > h && g[dq[t - 1]] >= g[k]) t--;
      dq[t++] = k;
      while (dq[h] > k + L - 1) h++;
      m[k] = g[dq[h]];
    }
    // moving average over m[i-L+1..i]
    var b = new Float32Array(n), acc = 0;
    for (var j = 0; j < n; j++) {
      acc += m[j]; if (j >= L) acc -= m[j - L];
      var cnt = Math.min(j + 1, L);
      b[j] = acc / cnt; // at the very start, average only what exists (keeps b <= g)
    }
    // release: gain may only recover slowly. min() with b keeps the brick-wall guarantee.
    var alpha = 1 - Math.exp(-1 / (relMs / 1000 * fs)), r = 1;
    for (var q = 0; q < n; q++) {
      r = r + (1 - r) * alpha;
      if (b[q] < r) r = b[q];
      b[q] = r;
    }
    return b;
  }

  function applyGain(chs, gainLin, curve) {
    return chs.map(function (x) {
      var y = new Float32Array(x.length);
      if (curve) for (var i = 0; i < x.length; i++) y[i] = x[i] * gainLin * curve[i];
      else for (var j = 0; j < x.length; j++) y[j] = x[j] * gainLin;
      return y;
    });
  }

  function db(x) { return 20 * Math.log10(x); }
  function lin(d) { return Math.pow(10, d / 20); }

  /* master(): chs = array of Float32Array, fs = sample rate.
     opts: targetLufs (e.g. -14), ceilingDb (e.g. -1). onProgress(text) optional.
     Returns { out, before, after, info } */
  function measure(chs, fs) {
    return { lufs: integratedLoudness(chs, fs), truePeak: truePeakDb(chs), samplePeak: samplePeakDb(chs), clipped: clippedSamples(chs) };
  }

  function master(chs, fs, opts, onProgress) {
    var step = onProgress || function () {};
    var target = opts.targetLufs, ceil = opts.ceilingDb;
    step('Measuring original');
    var before = measure(chs, fs);
    if (!isFinite(before.lufs)) throw new Error('This file looks silent.');
    var gainDb = target - before.lufs, usedLimiter = false, maxGr = 0, out, after, ceilTry = ceil - 0.1;
    var env = null;
    for (var iter = 0; iter < 6; iter++) {
      step('Processing, pass ' + (iter + 1));
      var g = lin(gainDb);
      var scaled = applyGain(chs, g);
      var peakDb = truePeakDb(scaled);
      if (peakDb <= ceil - 0.05) { out = scaled; usedLimiter = false; maxGr = 0; }
      else {
        usedLimiter = true;
        var e = peakEnvelope(scaled);
        var curve = limiterGain(e, lin(ceilTry), fs, 5, 120);
        var mn = 1; for (var i = 0; i < curve.length; i++) if (curve[i] < mn) mn = curve[i];
        maxGr = -db(mn);
        out = applyGain(scaled, 1, curve);
      }
      after = measure(out, fs);
      var err = target - after.lufs, tpOver = after.truePeak - ceil;
      if (tpOver > 0) { ceilTry -= tpOver + 0.05; if (Math.abs(err) < 0.1) continue; }
      if (Math.abs(err) <= 0.1 && tpOver <= 0) break;
      if (Math.abs(err) > 0.1) gainDb += err;
    }
    // Final check: how much of the target did we reach
    var reached = Math.abs(target - after.lufs) <= 0.2;
    return { out: out, before: before, after: after,
             info: { gainDb: after.lufs - before.lufs, usedLimiter: usedLimiter, maxGainReductionDb: maxGr, reachedTarget: reached } };
  }

  // ---- WAV writing ----
  function encodeWav(chs, fs, bits) {
    var nch = chs.length, n = chs[0].length, bytes = bits / 8, dataLen = n * nch * bytes;
    var buf = new ArrayBuffer(44 + dataLen), dv = new DataView(buf);
    function str(o, s) { for (var i = 0; i < s.length; i++) dv.setUint8(o + i, s.charCodeAt(i)); }
    str(0, 'RIFF'); dv.setUint32(4, 36 + dataLen, true); str(8, 'WAVE'); str(12, 'fmt ');
    dv.setUint32(16, 16, true); dv.setUint16(20, 1, true); dv.setUint16(22, nch, true);
    dv.setUint32(24, fs, true); dv.setUint32(28, fs * nch * bytes, true); dv.setUint16(32, nch * bytes, true);
    dv.setUint16(34, bits, true); str(36, 'data'); dv.setUint32(40, dataLen, true);
    var o = 44, max = Math.pow(2, bits - 1) - 1, prev = [0, 0];
    for (var i = 0; i < n; i++) for (var c = 0; c < nch; c++) {
      var v = chs[c][i] * max;
      if (bits === 16) { var r = Math.random(); v += (r - prev[c]); prev[c] = r; } // TPDF dither (high-passed)
      v = Math.round(v); if (v > max) v = max; if (v < -max - 1) v = -max - 1;
      if (bits === 16) dv.setInt16(o, v, true);
      else { dv.setUint8(o, v & 255); dv.setUint8(o + 1, (v >> 8) & 255); dv.setUint8(o + 2, (v >> 16) & 255); }
      o += bytes;
    }
    return buf;
  }
  function wavSampleRate(arrayBuffer) {
    try {
      var dv = new DataView(arrayBuffer);
      if (dv.getUint32(0) === 0x52494646 && dv.getUint32(8) === 0x57415645) return dv.getUint32(24, true);
    } catch (e) {}
    return null;
  }

  var api = { integratedLoudness: integratedLoudness, truePeakDb: truePeakDb, samplePeakDb: samplePeakDb,
              clippedSamples: clippedSamples, measure: measure, master: master, encodeWav: encodeWav, wavSampleRate: wavSampleRate,
              _limiterGain: limiterGain, _peakEnvelope: peakEnvelope, _applyGain: applyGain };
  if (typeof module !== 'undefined' && module.exports) module.exports = api; else root.AutoMaster = api;
})(typeof self !== 'undefined' ? self : this);
