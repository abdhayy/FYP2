import React, { useMemo } from "react";
import {
  Area, Bar, BarChart, CartesianGrid, Cell, ComposedChart, Legend, Line, LineChart, ReferenceArea, ReferenceLine,
  ResponsiveContainer, Tooltip, XAxis, YAxis,
} from "recharts";
import { Panel } from "./Window1.jsx";
import { GREEN, ORANGE, RED, fmt, signed } from "./common.js";

const NO_ANIM = { isAnimationActive: false };

export default function Window2({ s, w }) {
  const E = s.entropy;
  const d = w.w2;
  const zc = E.zClear;
  const faulty = d.verdict === "FAULTY";
  const clear = Object.fromEntries(E.keys.map((k) => [k, Math.abs(d.z[k]) >= zc]));
  const nClear = E.keys.filter((k) => clear[k]).length;
  const nRel = E.keys.filter((k) => clear[k] && E.reliable.includes(k)).length;
  const short = (k) => E.names[k].replace(" Entropy", "");

  // the gauge range only grows over the windows already played, like the Python version
  const half = useMemo(() => {
    const h = {};
    for (const k of E.keys) {
      h[k] = 8;
      for (const sw of s.windows) if (sw.index <= w.index) h[k] = Math.max(h[k], Math.abs(sw.z[k] ?? 0) + 1.5);
    }
    return h;
  }, [s, w.index, E.keys]);

  const barRows = [...E.keys]
    .sort((a, b) => Math.abs(d.z[b]) - Math.abs(d.z[a]))
    .map((k) => ({ k: short(k), v: Math.abs(d.z[k]), c: clear[k] ? (E.reliable.includes(k) ? RED : ORANGE) : "#d3d3d3" }));

  const trendRows = s.windows.filter((sw) => sw.index <= w.index).map((sw) => ({ n: sw.index, ...sw.z }));
  const allZ = trendRows.flatMap((r) => E.keys.map((k) => r[k] ?? 0));
  const zLo = Math.min(-6, Math.min(...allZ) - 1), zHi = Math.max(6, Math.max(...allZ) + 1);

  const mseRows = d.mse.map((v, i) => ({
    sc: i + 1, v, mean: E.mseMean[i], band: [E.mseMean[i] - 3 * E.mseStd[i], E.mseMean[i] + 3 * E.mseStd[i]],
  }));
  const verdictColors = ["#1f77b4", "#ff7f0e"];

  return (
    <div className="window">
      <div className="banner-wrap">
        <div className="win-title">
          WINDOW 2 — ENTROPY: IS THERE A FAULT, AND WHICH ENTROPY SHOWS IT · #{w.index} · window {fmt(w.tStart, 0)}–{fmt(w.tEnd, 0)} s · {s.label}
          <span className="muted"> [computed in {fmt(d.elapsed, 1)} s, entropies at {E.analysisFs} Hz]</span>
        </div>
        <div className="banner" style={{ color: faulty ? RED : GREEN }}>ENTROPY VERDICT: {d.verdict}</div>
        <div className="sub">
          {nClear} of 7 entropies are clearly outside the healthy range now ({nRel} of them reliable) | rule: FAULTY when{" "}
          {E.verdictKeys.map(short).join(" AND ")} are both below their cut-off
        </div>
      </div>

      <div className="grid-e">
        <div className="panel">
          <div className="panel-title bold navy">
            EACH ENTROPY vs A HEALTHY BEARING (green bar = normal healthy range, red dashed = cut-off, dot = this window)
          </div>
          {E.keys.map((k) => (
            <Gauge key={k} k={k} name={E.names[k]} value={d.vals[k]} z={d.z[k]} mean={E.mean[k]} std={E.std[k]}
              half={half[k]} zc={zc} clear={clear[k]} reliable={E.reliable.includes(k)} test={E.test[k]} />
          ))}
        </div>

        <div className="col">
          <Panel title="This window: how far from healthy — red = reliable, orange = not reliable" bold>
            <ResponsiveContainer width="100%" height={210}>
              <BarChart data={barRows} layout="vertical" margin={{ top: 4, right: 12, bottom: 18, left: 8 }}>
                <CartesianGrid strokeOpacity={0.3} horizontal={false} />
                <XAxis type="number" domain={[0, (m) => Math.max(4, Math.ceil(m))]} label={{ value: "how far from healthy (std. dev.)", position: "insideBottom", offset: -8, fontSize: 11 }} />
                <YAxis type="category" dataKey="k" width={86} tick={{ fontSize: 11 }} />
                <Tooltip formatter={(v) => fmt(v, 2)} />
                <ReferenceLine x={zc} stroke={RED} strokeDasharray="5 3" />
                <Bar dataKey="v" stroke="#000" strokeWidth={0.5} {...NO_ANIM}>
                  {barRows.map((r) => <Cell key={r.k} fill={r.c} />)}
                </Bar>
              </BarChart>
            </ResponsiveContainer>
          </Panel>

          <Panel title="Trend (z-score per 10 s window)" bold>
            <ResponsiveContainer width="100%" height={200}>
              <LineChart data={trendRows} margin={{ top: 4, right: 12, bottom: 18, left: 0 }}>
                <CartesianGrid strokeOpacity={0.3} />
                <XAxis dataKey="n" type="number" domain={[0.5, s.windows.length + 0.5]} ticks={[1, ...Array.from({ length: Math.floor(s.windows.length / 5) }, (_, i) => (i + 1) * 5)]}
                  label={{ value: "10 s window number", position: "insideBottom", offset: -8, fontSize: 11 }} />
                <YAxis width={36} domain={[Math.floor(zLo), Math.ceil(zHi)]} />
                <Tooltip formatter={(v) => fmt(v, 2)} labelFormatter={(n) => `window ${n}`} />
                <ReferenceArea y1={-zc} y2={zc} fill={GREEN} fillOpacity={0.08} />
                <ReferenceLine y={-zc} stroke={RED} strokeDasharray="5 3" />
                <ReferenceLine y={zc} stroke={RED} strokeDasharray="5 3" />
                {E.keys.filter((k) => !E.verdictKeys.includes(k)).map((k) => (
                  <Line key={k} dataKey={k} name={short(k)} stroke="#999" strokeOpacity={0.6} strokeWidth={0.8} dot={false} legendType="none" {...NO_ANIM} />
                ))}
                {E.verdictKeys.map((k, i) => (
                  <Line key={k} dataKey={k} name={short(k)} stroke={verdictColors[i]} strokeWidth={1.5} dot={{ r: 2 }} {...NO_ANIM} />
                ))}
                <Legend verticalAlign="top" height={20} wrapperStyle={{ fontSize: 11 }} />
              </LineChart>
            </ResponsiveContainer>
          </Panel>

          <Panel title="Multiscale Sample Entropy" bold>
            <ResponsiveContainer width="100%" height={200}>
              <ComposedChart data={mseRows} margin={{ top: 4, right: 12, bottom: 18, left: 0 }}>
                <CartesianGrid strokeOpacity={0.3} />
                <XAxis dataKey="sc" type="number" domain={[1, mseRows.length]} ticks={mseRows.map((r) => r.sc)} label={{ value: "scale factor", position: "insideBottom", offset: -8, fontSize: 11 }} />
                <YAxis width={36} tickFormatter={(v) => v.toFixed(1)} />
                <Tooltip formatter={(v) => (Array.isArray(v) ? `${fmt(v[0], 3)} – ${fmt(v[1], 3)}` : fmt(v, 3))} />
                <Area dataKey="band" name="healthy range" fill={GREEN} fillOpacity={0.22} stroke="none" {...NO_ANIM} />
                <Line dataKey="mean" name="healthy mean" stroke={GREEN} dot={false} strokeWidth={1} legendType="none" {...NO_ANIM} />
                <Line dataKey="v" name="this window" stroke={faulty ? RED : "#1f77b4"} strokeWidth={1.6} dot={{ r: 3 }} {...NO_ANIM} />
                <Legend verticalAlign="top" height={20} wrapperStyle={{ fontSize: 11 }} />
              </ComposedChart>
            </ResponsiveContainer>
          </Panel>
        </div>
      </div>

      <p className="howto">
        HOW TO READ: entropy = how irregular the vibration is. A faulty bearing repeats impacts, so the signal becomes more
        regular and Dispersion / Fuzzy entropy DROP below the healthy range. Healthy range: {E.source}; cut-off = healthy
        mean ± {zc} × std. "tested" = how many windows of your four recordings each entropy flagged (fixed note from the
        10 Nov test). Entropy says THAT a fault exists, not WHICH one (windows 1 and 3).
      </p>
    </div>
  );
}

// One horizontal gauge: pink = outside the healthy range, green bar = mean ± zc·std, dot = this window.
function Gauge({ name, value, z, mean, std, half, zc, clear, reliable, test }) {
  const lo = mean - half * std, hi = mean + half * std;
  const pct = (v) => `${Math.min(100, Math.max(0, ((v - lo) / (hi - lo)) * 100))}%`;
  const cLo = mean - zc * std, cHi = mean + zc * std;
  const tag = clear && reliable ? "RELIABLE fault signal" : clear ? "changed here, but NOT reliable" : "no clear change";
  const tc = clear && reliable ? RED : clear ? ORANGE : "#666";
  const ticks = niceTicks(lo, hi, 6);
  return (
    <div className="gauge">
      <div className="g-name">{name}</div>
      <div className="g-track-wrap">
        <div className="g-track">
          <div className="g-out" style={{ left: 0, width: pct(cLo) }} />
          <div className="g-out" style={{ left: pct(cHi), right: 0 }} />
          <div className="g-ok" style={{ left: pct(cLo), width: `calc(${pct(cHi)} - ${pct(cLo)})` }} />
          <div className="g-cut" style={{ left: pct(cLo) }} />
          <div className="g-cut" style={{ left: pct(cHi) }} />
          <div className="g-dot" style={{ left: pct(value), background: clear ? RED : "#2ca02c" }} title={`${fmt(value, 4)}`} />
        </div>
        <div className="g-axis">
          {ticks.map((t) => (
            <span key={t} style={{ left: pct(t) }}>{+t.toFixed(4)}</span>
          ))}
        </div>
      </div>
      <div className="g-text" style={{ color: tc, fontWeight: clear && reliable ? 700 : 400 }}>
        {fmt(value, 3)} &nbsp; z = {signed(z, 1)} ({z < 0 ? "lower" : "higher"})
        <br />
        {tag}
        <br />
        tested: faulty {test[0]}/60, healthy {test[1]}/30
      </div>
    </div>
  );
}

function niceTicks(lo, hi, n) {
  const span = hi - lo;
  const raw = span / n;
  const mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => s >= raw) || raw;
  const out = [];
  for (let t = Math.ceil(lo / step) * step; t <= hi + 1e-12; t += step) out.push(Number(t.toPrecision(10)));
  return out;
}
