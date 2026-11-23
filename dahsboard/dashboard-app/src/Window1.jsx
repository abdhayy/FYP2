import React, { useMemo } from "react";
import {
  Area, Bar, BarChart, Cell, ComposedChart, LabelList, Legend, Line, LineChart, ReferenceDot, ReferenceLine,
  ResponsiveContainer, Tooltip, XAxis, YAxis, CartesianGrid,
} from "recharts";
import { FAULTS, FAULT_COLORS, GREEN, RED, dashOf, envMarkers, fftMarkers, fmt, zipRows } from "./common.js";

const NO_ANIM = { isAnimationActive: false, dot: false };

export default function Window1({ s, w }) {
  const d = w.w1;
  const faulty = d.verdict !== "HEALTHY";
  const thr = s.faultDb;

  const timeRows = useMemo(() => zipRows({ t: d.time.t, raw: d.time.raw, smooth: d.time.smooth }), [d]);
  const fftRows = useMemo(() => zipRows({ f: d.fft.f, a: d.fft.a }), [d]);
  const compRows = useMemo(() => {
    const cols = { t: d.compT };
    FAULTS.forEach((k, i) => {
      const peak = Math.max(...d.components[k].y.map(Math.abs)) || 1;
      cols[k] = d.components[k].y.map((v) => 3 - i + (0.42 * v) / peak);
    });
    return zipRows(cols);
  }, [d]);
  const sumRows = useMemo(() => zipRows({ t: d.compT, env: d.sum.env, sum: d.sum.sum }), [d]);
  const envRows = useMemo(() => {
    const ref = s.refSpectrum && s.refSpectrum.a.length === d.env.f.length ? s.refSpectrum.a : null;
    return d.env.f.map((f, i) => {
      const a = Math.max(d.env.a[i], 1e-9);
      const t = d.env.thr[i];
      return {
        f, a, floor: d.env.floor[i], thr: t, ref: ref ? Math.max(ref[i], 1e-9) : undefined,
        ev: a > t ? [t, a] : [t, t],
      };
    });
  }, [d, s]);
  const envDomain = useMemo(() => {
    // like log_limits() in the Python script, but robust: ignore the roll-off near 450 Hz
    const v = d.env.a.filter((x, i) => d.env.f[i] > 3 && d.env.f[i] < 400 && x > 0).sort((a, b) => a - b);
    const lo = v[Math.floor(v.length * 0.02)] / 4, hi = v[v.length - 1] * 3;
    return [Math.pow(10, Math.floor(Math.log10(lo))), Math.pow(10, Math.ceil(Math.log10(hi)))];
  }, [d]);
  const evRows = FAULTS.map((k) => ({ k, v: d.evidence[k], known: s.known[k] }));

  return (
    <div className="window">
      <div className="banner-wrap">
        <div className="banner" style={{ color: faulty ? RED : GREEN }}>CONDITION: {d.verdict}</div>
        <div className="sub">
          evidence (dB above the noise floor; fault threshold {fmt(thr, 2)} dB):{" "}
          {FAULTS.map((k) => (
            <span key={k} className="ev-chip" style={{ color: FAULT_COLORS[k] }}>
              {k} {fmt(d.evidence[k], 0)}
              {s.known[k] ? "" : "*"}
            </span>
          ))}
          {d.others.length > 0 && <span> | also present: {d.others.join(", ")}</span>}
        </div>
      </div>

      <div className="grid2">
        <Panel title={<>Time Domain Waveform (Fs = {s.fs} Hz) — last {s.groupSeconds} s{" "}
          <span className={"status " + (d.fft.detected ? "ok" : "")}>{d.fft.detected ? "VIBRATION DETECTED" : "IDLE – no significant vibration"}</span></>}>
          <ResponsiveContainer width="100%" height={230}>
            <LineChart data={timeRows} margin={{ top: 8, right: 12, bottom: 18, left: 8 }}>
              <CartesianGrid strokeOpacity={0.4} />
              <XAxis dataKey="t" type="number" domain={[0, s.groupSeconds]} tickCount={6} label={{ value: "Time (s)", position: "insideBottom", offset: -8 }} />
              <YAxis width={72} tickFormatter={(v) => v.toFixed(2)} label={{ value: "Amplitude (g)", angle: -90, position: "insideLeft", offset: 0, fontSize: 11 }} />
              <Tooltip formatter={(v) => fmt(v, 4)} labelFormatter={(t) => `t = ${fmt(t, 3)} s`} />
              <Legend verticalAlign="top" height={22} />
              <Line dataKey="raw" name="Raw" stroke="#b0c4de" strokeWidth={0.8} {...NO_ANIM} />
              <Line dataKey="smooth" name="Smoothed (150 Hz LP)" stroke="#0000ff" strokeWidth={1.1} {...NO_ANIM} />
            </LineChart>
          </ResponsiveContainer>
        </Panel>

        <Panel title="FFT Spectrum (3 Harmonic Modes) — dotted lines = BPFO, BPFI, 2xBSF, FTF">
          <ResponsiveContainer width="100%" height={250}>
            <LineChart data={fftRows} margin={{ top: 18, right: 16, bottom: 18, left: 8 }}>
              <CartesianGrid strokeOpacity={0.4} />
              <XAxis dataKey="f" type="number" domain={[0, 300]} ticks={[0, 50, 100, 150, 200, 250, 300]} label={{ value: "Frequency (Hz)", position: "insideBottom", offset: -8 }} />
              <YAxis width={72} tickFormatter={(v) => v.toFixed(3)} label={{ value: "Amplitude (g)", angle: -90, position: "insideLeft", offset: 0, fontSize: 11 }} />
              <Tooltip formatter={(v) => fmt(v, 5)} labelFormatter={(f) => `${fmt(f, 1)} Hz`} />
              {fftMarkers(s.bearing, s.known).map((m) => (
                <ReferenceLine key={m.name} x={m.f} stroke={FAULT_COLORS[m.key]} strokeDasharray="2 3" strokeWidth={1.3}
                  label={{ value: `${m.name} ${fmt(m.f, 1)}`, angle: -90, position: "insideBottomLeft", fill: FAULT_COLORS[m.key], fontSize: 10 }} />
              ))}
              <Line dataKey="a" stroke="#e00000" strokeWidth={1.2} {...NO_ANIM} />
              {d.fft.modes.map((m) => (
                <ReferenceDot key={m.k} x={m.f} y={m.a} r={6} fill="#ffff00" stroke="#000"
                  label={{ value: `Mode ${m.k}: ${fmt(m.f, 1)} Hz`, position: "right", fontSize: 11 }} />
              ))}
            </LineChart>
          </ResponsiveContainer>
        </Panel>

        <Panel title="TIME DOMAIN — the four fault components extracted from YOUR signal (envelope domain, each drawn to its own peak)">
          <ResponsiveContainer width="100%" height={260}>
            <LineChart data={compRows} margin={{ top: 8, right: 12, bottom: 18, left: 8 }}>
              <CartesianGrid strokeOpacity={0.4} />
              <XAxis dataKey="t" type="number" domain={[0, "dataMax"]} tickFormatter={(v) => v.toFixed(1)} label={{ value: "Time (s) — last 1 s", position: "insideBottom", offset: -8 }} />
              <YAxis width={92} domain={[-0.6, 3.6]} ticks={[0, 1, 2, 3]}
                tick={({ x, y, payload }) => {
                  const k = FAULTS[3 - payload.value];
                  return (
                    <g transform={`translate(${x},${y})`}>
                      <text x={-4} y={-2} textAnchor="end" fontSize={11} fill={FAULT_COLORS[k]} fontWeight="bold">{k}{s.known[k] ? "" : "*"}</text>
                      <text x={-4} y={11} textAnchor="end" fontSize={10} fill="#555">peak {fmt(d.components[k].peak * 1000, 1)} mg</text>
                    </g>
                  );
                }} />
              {FAULTS.map((k) => (
                <Line key={k} dataKey={k} stroke={FAULT_COLORS[k]} strokeWidth={1} {...NO_ANIM} />
              ))}
            </LineChart>
          </ResponsiveContainer>
        </Panel>

        <Panel title="WHAT DEFINES THE FAULT — envelope FFT: lines rising above the red threshold" bold>
          <ResponsiveContainer width="100%" height={300}>
            <ComposedChart data={envRows} margin={{ top: 18, right: 16, bottom: 18, left: 8 }}>
              <CartesianGrid strokeOpacity={0.4} />
              <XAxis dataKey="f" type="number" domain={[0, 450]} ticks={[0, 50, 100, 150, 200, 250, 300, 350, 400, 450]} label={{ value: "Frequency (Hz)", position: "insideBottom", offset: -8 }} />
              <YAxis width={62} scale="log" domain={envDomain} allowDataOverflow tickFormatter={(v) => v.toExponential(0)}
                label={{ value: "Envelope (g, log)", angle: -90, position: "insideLeft", offset: 0, fontSize: 11 }} />
              <Tooltip formatter={(v) => (Array.isArray(v) ? null : Number(v).toExponential(2))} labelFormatter={(f) => `${fmt(f, 2)} Hz`} />
              <Legend verticalAlign="bottom" height={34} wrapperStyle={{ fontSize: 11, bottom: -6 }} />
              {envMarkers(s.bearing).map((m) => (
                <ReferenceLine key={m.name} x={m.f} stroke={FAULT_COLORS[m.fault]} strokeDasharray={dashOf(m.ls)} strokeOpacity={0.85}
                  label={m.ls === "dot" ? undefined : { value: `${m.name} ${fmt(m.f, 1)}`, angle: -90, position: m.ls === "solid" ? "insideTopRight" : "insideRight", fill: FAULT_COLORS[m.fault], fontSize: 9 }} />
              ))}
              <ReferenceLine x={s.bearing.fr} stroke="#888" strokeDasharray="6 2 2 2" />
              <Area dataKey="ev" name="evidence (above threshold)" stroke="none" fill={RED} fillOpacity={0.35} isAnimationActive={false} legendType="square" />
              <Line dataKey="a" name="envelope spectrum of YOUR signal" stroke="#000" strokeWidth={0.9} {...NO_ANIM} />
              {s.refSpectrum && <Line dataKey="ref" name="Healthy recording (reference)" stroke="#1a991a" strokeWidth={0.9} strokeOpacity={0.85} {...NO_ANIM} />}
              <Line dataKey="floor" name="noise floor (local median)" stroke="#696969" strokeDasharray="5 3" strokeWidth={1.4} {...NO_ANIM} />
              <Line dataKey="thr" name={`fault threshold = floor + ${fmt(thr, 0)} dB`} stroke={RED} strokeWidth={1.2} {...NO_ANIM} />
              {d.env.marks.filter((m) => m.above).map((m) => (
                <ReferenceDot key={m.name} x={m.f} y={Math.max(m.a, envDomain[0])} r={5} fill={FAULT_COLORS[m.fault]} stroke="#000"
                  label={{ value: `${m.name}${s.known[m.fault] ? "" : "*"} ${fmt(m.f, 1)} Hz +${fmt(m.db, 0)} dB`, position: "top", fill: FAULT_COLORS[m.fault], fontSize: 10, fontWeight: "bold" }} />
              ))}
            </ComposedChart>
          </ResponsiveContainer>
        </Panel>

        <Panel title="TIME DOMAIN — the four components COMBINED (sum)" bold>
          <ResponsiveContainer width="100%" height={230}>
            <LineChart data={sumRows} margin={{ top: 8, right: 12, bottom: 18, left: 8 }}>
              <CartesianGrid strokeOpacity={0.4} />
              <XAxis dataKey="t" type="number" domain={[0, "dataMax"]} tickFormatter={(v) => v.toFixed(1)} label={{ value: "Time (s) — last 1 s", position: "insideBottom", offset: -8 }} />
              <YAxis width={72} tickFormatter={(v) => v.toFixed(3)} label={{ value: "Envelope (g)", angle: -90, position: "insideLeft", offset: 0, fontSize: 11 }} />
              <Tooltip formatter={(v) => fmt(v, 5)} labelFormatter={(t) => `t = ${fmt(t, 3)} s`} />
              <Legend verticalAlign="top" height={22} />
              <Line dataKey="env" name="measured envelope" stroke="#c8c8c8" strokeWidth={1.4} {...NO_ANIM} />
              <Line dataKey="sum" name="COMBINED = BPFO + BPFI + BSF + FTF" stroke="#000" strokeWidth={1} {...NO_ANIM} />
            </LineChart>
          </ResponsiveContainer>
        </Panel>

        <div className="grid-inner">
          <Panel title="Evidence per fault (dB above noise floor)">
            <ResponsiveContainer width="100%" height={230}>
              <BarChart data={evRows} margin={{ top: 18, right: 8, bottom: 4, left: 0 }}>
                <CartesianGrid strokeOpacity={0.3} vertical={false} />
                <XAxis dataKey="k" />
                <YAxis width={34} domain={[0, (max) => Math.max(30, Math.ceil(max * 1.25 / 10) * 10)]} tickCount={4} allowDecimals={false} />
                <Tooltip formatter={(v) => `${fmt(v, 1)} dB`} />
                <ReferenceLine y={thr} stroke={RED} strokeWidth={1.4} label={{ value: `threshold ${fmt(thr, 1)} dB`, position: "insideTopRight", fill: RED, fontSize: 10 }} />
                <Bar dataKey="v" isAnimationActive={false} stroke="#000" strokeWidth={0.6}>
                  {evRows.map((r) => <Cell key={r.k} fill={FAULT_COLORS[r.k]} fillOpacity={r.known ? 1 : 0.45} />)}
                  <LabelList dataKey="v" position="top" formatter={(v) => fmt(v, 0)} style={{ fontWeight: "bold", fontSize: 12 }} />
                </Bar>
              </BarChart>
            </ResponsiveContainer>
          </Panel>
          <RuleBox s={s} d={d} />
        </div>
      </div>
    </div>
  );
}

function RuleBox({ s, d }) {
  const b = s.bearing;
  const thr = s.faultDb;
  const faulty = d.verdict !== "HEALTHY";
  const lines = [
    "HOW THE VERDICT IS DECIDED",
    `1 band-pass ${s.band[0]}-${s.band[1]} Hz -> envelope -> spectrum`,
    "2 each fault has fixed LINE frequencies:",
    `  FTF  = cage  ${fmt(b.FTF)} Hz and x2 x3 x4`,
    `  BSF  = 2xBSF ${fmt(2 * b.BSF)} Hz (+/-FTF)`,
    `  BPFO = ${fmt(b.BPFO)}${s.known.BPFO ? "" : "*"} Hz   BPFI = ${fmt(b.BPFI)}${s.known.BPFI ? "" : "*"} Hz`,
    "3 evidence = mean dB of the 2 strongest",
    "  lines above the local noise floor",
    `4 FAULT if evidence >= ${fmt(thr, 1)} dB`,
    "  (Healthy recordings never exceed it)",
    "-".repeat(40),
    ...(faulty
      ? [
          `NOW: ${d.top} lines are ${fmt(d.topDb, 0)} dB high:`,
          ...d.topLines.map((l) => `  ${l.name} ${fmt(l.f)} Hz  +${fmt(l.db, 0)} dB ${l.db >= thr ? "<- above" : ""}`),
          `=> ${d.top} FAULTY (${fmt(d.topDb, 0)} >= ${fmt(thr, 1)} dB)`,
        ]
      : [`NOW: strongest = ${d.top} ${fmt(d.topDb, 0)} dB`, `  < ${fmt(thr, 1)} dB threshold`, "=> HEALTHY (no fault line rises out of", "   the noise floor)"]),
  ];
  return (
    <div className="panel">
      <pre className="rulebox" style={{ color: faulty ? RED : "#007300" }}>{lines.join("\n")}</pre>
    </div>
  );
}

export function Panel({ title, children, bold }) {
  return (
    <div className="panel">
      <div className={"panel-title" + (bold ? " bold" : "")}>{title}</div>
      {children}
    </div>
  );
}
