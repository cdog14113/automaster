# Auto Master

Browser-only tool that measures a song's loudness and true peak, then sets it to a target (default -14 LUFS, -1 dBTP) with plain gain, adding a lookahead limiter only if peaks need it. Audio never leaves the device.

- `dsp.js`: BS.1770 loudness, 4x true-peak, limiter, WAV writer
- `index.html`: the UI
- `test.js`, `test2.js`: Node checks (results were compared with ffmpeg's ebur128 meter)
