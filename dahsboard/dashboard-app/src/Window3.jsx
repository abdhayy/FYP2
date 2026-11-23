import React, { useMemo } from "react";
import {
  CartesianGrid, Legend, Line, LineChart, ReferenceLine, ResponsiveContainer, Tooltip, XAxis, YAxis,
} from "recharts";
import { Panel } from "./Window1.jsx";
import { FAULTS, FAULT_COLORS, GREEN, RED, fmt, zipRows } from "./common.js";

const NO_ANIM = { isAnimationActive: false, dot: false };

export default function Window3({ s, w }) {
  const d = w.w3;
  const faulty = d.verdict !== "HEALTHY";
  const combRows = useMemo(() => zipRows({ t: d.combined.t, sig: d.combined.sig, sum: d.combined.modesSum }), [d]);
  const specRows = useMemo(() => zipRows({ f: d.spectrum.f, a: d.spectrum.a }), [d]);

  return (
    <div className="window">
      <div className="banner-wrap">
        <div className="win-title">
          WINDOW 3 — VMD OF THE COMBINED TIME-DOMAIN SIGNAL · #{w.index} · window {fmt(w.tStart, 0)}–{fmt(w.tEnd, 0)} s · K = {d.K} modes · {s.label}
          <span className="muted"> [computed in {fmt(d.elapsed, 1)} s, {d.nIter} iterations]</span>
        </div>
        <div className="banner" style={{ color: faulty ? RED : GREEN }}>CONDITION: {d.verdict}</div>
      </div>

      <div className="grid-v">
        <div className="col">
          <Panel title={`Time domain: BPFO + BPFI + BSF + FTF combined → VMD   (sum of the ${d.K} VMD modes, residual ${fmt(d.resid * 100, 1)} %)`} bold>
            <ResponsiveContainer width="100%" height={170}>
              <LineChart data={combRows} margin={{ top: 4, right: 12, bottom: 4, left: 8 }}>
                <CartesianGrid strokeOpacity={0.3} />
                <XAxis dataKey="t" type="number" domain={[0, "dataMax"]} tickFormatter={(v) => v.toFixed(2)} />
                <YAxis width={58} tickFormatter={(v) => v.toFixed(3)} />
                <Tooltip formatter={(v) => fmt(v, 5)} labelFormatter={(t) => `t = ${fmt(t, 3)} s`} />
                <Legend verticalAlign="top" height={20} wrapperStyle={{ fontSize: 11 }} />
                <Line dataKey="sig" name="COMBINED signal (input of VMD)" stroke="#000" strokeWidth={1.1} {...NO_ANIM} />
                <Line dataKey="sum" name={`sum of the ${d.K} VMD modes`} stroke="#ff7f0e" strokeDasharray="4 3" strokeWidth={1.1} {...NO_ANIM} />
              </LineChart>
            </ResponsiveContainer>
          </Panel>
          {FAULTS.map((k) => <FaultPanel key={k} k={k} f={d.faults[k]} t={d.combined.t} last={k === "FTF"} />)}
        </div>

        <div className="col">
          <Panel title="Spectrum of the combined signal — each VMD mode sits on one fault line (colour = fault)" bold>
            <ResponsiveContainer width="100%" height={300}>
              <LineChart data={specRows} margin={{ top: 70, right: 16, bottom: 18, left: 8 }}>
                <CartesianGrid strokeOpacity={0.3} />
                <XAxis dataKey="f" type="number" domain={[0, 330]} ticks={[0, 50, 100, 150, 200, 250, 300]} label={{ value: "Frequency (Hz)", position: "insideBottom", offset: -8 }} />
                <YAxis width={44} tickFormatter={(v) => v.toFixed(1)} label={{ value: "mg", angle: -90, position: "insideLeft" }} />
                <Tooltip formatter={(v) => `${fmt(v, 3)} mg`} labelFormatter={(f) => `${fmt(f, 2)} Hz`} />
                {d.modes.map((m) => (
                  <ReferenceLine key={m.i} x={m.center} stroke={FAULT_COLORS[m.fault] || "#999"} strokeDasharray="2 3"
                    label={{ value: `IMF${m.i} ${m.line || "?"}`, angle: -90, position: "insideTopLeft", fill: FAULT_COLORS[m.fault] || "#999", fontSize: 9 }} />
                ))}
                <Line dataKey="a" stroke="#000" strokeWidth={1} {...NO_ANIM} />
              </LineChart>
            </ResponsiveContainer>
          </Panel>

          <div className="panel">
            <div className="panel-title bold">VMD modes and the fault line each one sits on</div>
            <table className="tbl">
              <thead>
                <tr><th>IMF</th><th>centre Hz</th><th>fault : line</th><th>RMS mg</th></tr>
              </thead>
              <tbody>
                {d.modes.map((m) => (
                  <tr key={m.i}>
                    <td>{m.i}</td>
                    <td>{fmt(m.center, 1)}</td>
                    <td style={{ color: FAULT_COLORS[m.fault] || "#999", fontWeight: 600 }}>{m.fault ? `${m.fault} : ${m.line}` : "– (unlabelled)"}</td>
                    <td>{fmt(m.rms, 2)}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>

          <div className="panel">
            <div className="panel-title bold">FAULT FEATURES</div>
            <table className="tbl">
              <thead>
                <tr><th>fault</th><th>evidence dB</th><th>RMS mg</th><th>align r</th><th>share %</th><th>IMFs</th></tr>
              </thead>
              <tbody>
                {FAULTS.map((k) => {
                  const f = d.faults[k];
                  return (
                    <tr key={k}>
                      <td style={{ color: FAULT_COLORS[k], fontWeight: 700 }}>{k}</td>
                      <td style={{ fontWeight: f.evidence >= s.faultDb ? 700 : 400, color: f.evidence >= s.faultDb ? RED : undefined }}>{fmt(f.evidence, 1)}</td>
                      <td>{fmt(f.rms, 2)}</td>
                      <td>{fmt(f.r, 3)}</td>
                      <td>{fmt(f.share * 100, 0)}</td>
                      <td>{f.imfs.join(", ")}</td>
                    </tr>
                  );
                })}
              </tbody>
            </table>
          </div>
        </div>
      </div>
      <p className="howto">
        VMD recovering the four parts shows the decomposition works (r ≈ 0.99 also on Healthy data). The proof that a fault
        exists is the evidence in dB against the noise floor (window 1); entropy is in window 2.
      </p>
    </div>
  );
}

function FaultPanel({ k, f, t, last }) {
  const rows = useMemo(() => zipRows({ t, orig: f.orig, rec: f.rec }), [t, f]);
  return (
    <div className="panel">
      <div className="panel-title bold" style={{ color: FAULT_COLORS[k] }}>
        {k} = IMF {f.imfs.join(", ")} &nbsp; r = {fmt(f.r, 3)} &nbsp; RMS {fmt(f.rms, 2)} mg &nbsp; share {fmt(f.share * 100, 0)} % &nbsp; evidence {fmt(f.evidence, 0)} dB
      </div>
      <ResponsiveContainer width="100%" height={last ? 140 : 120}>
        <LineChart data={rows} margin={{ top: 2, right: 12, bottom: last ? 16 : 2, left: 8 }}>
          <CartesianGrid strokeOpacity={0.3} />
          <XAxis dataKey="t" type="number" domain={[0, "dataMax"]} tickFormatter={(v) => v.toFixed(2)} hide={!last}
            label={last ? { value: "Time (s) — last 2 s of the window", position: "insideBottom", offset: -8, fontSize: 11 } : undefined} />
          <YAxis width={58} tickFormatter={(v) => v.toFixed(4)} />
          <Tooltip formatter={(v) => fmt(v, 5)} labelFormatter={(x) => `t = ${fmt(x, 3)} s`} />
          <Line dataKey="orig" name="original component (before VMD)" stroke="#bbb" strokeWidth={3} {...NO_ANIM} />
          <Line dataKey="rec" name="separated again by VMD" stroke={FAULT_COLORS[k]} strokeWidth={1} {...NO_ANIM} />
        </LineChart>
      </ResponsiveContainer>
    </div>
  );
}
