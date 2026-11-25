import React, { useEffect, useMemo, useState } from "react";
import { CartesianGrid, ReferenceLine, ResponsiveContainer, Scatter, ScatterChart, Tooltip, XAxis, YAxis, ZAxis, Legend } from "recharts";
import { GREEN, RED, fmt } from "./common.js";

const CLASSES = ["normal", "cage_fault", "outer_race", "ball_fault"];
const NAMES = { normal: "Normal", cage_fault: "Cage (FTF)", outer_race: "Outer race (BPFO)", ball_fault: "Ball (BSF)" };
const COLORS = { normal: "#2ca02c", cage_fault: "#d62728", outer_race: "#ff7f0e", ball_fault: "#9467bd" };
const pct = (v) => (v === undefined || v === null ? "–" : `${Math.round(v * 100)} %`);

// Tab with the results of mafaulda_analysis.py (public MaFaulDa dataset, UFRJ).
export default function MaFaulda() {
  const [d, setD] = useState(null);
  const [missing, setMissing] = useState(false);
  const [zKey, setZKey] = useState("cage_fault");

  useEffect(() => {
    fetch("./data/mafaulda.json")
      .then((r) => (r.ok ? r.json() : Promise.reject()))
      .then(setD)
      .catch(() => setMissing(true));
  }, []);

  const scatter = useMemo(() => {
    if (!d) return {};
    const out = {};
    for (const c of CLASSES) out[c] = d.files_detail.filter((f) => f.cls === c && f.split === "test" && f.z[zKey] !== null).map((f) => ({ fr: f.fr, z: f.z[zKey], file: f.file }));
    return out;
  }, [d, zKey]);

  if (missing)
    return (
      <div className="window">
        <div className="panel">
          <div className="panel-title bold">MaFaulDa — no results yet</div>
          <p>Run the analysis on your laptop (where the dataset is), from the <code>dahsboard</code> folder:</p>
          <pre className="rulebox">python mafaulda_analysis.py --root "C:\Users\User\Desktop\python fyp"</pre>
          <p className="muted small">It writes <code>dashboard-app/public/data/mafaulda.json</code>, which this tab reads.</p>
        </div>
      </div>
    );
  if (!d) return <div className="card loading">Loading MaFaulDa results…</div>;

  const E = d.envelope_choice, N = d.entropy_choice;
  return (
    <div className="window">
      <div className="banner-wrap">
        <div className="win-title">
          MaFaulDa — the FYP method tested on a public dataset (UFRJ Machinery Fault Database, 50 kHz, 8-ball bearings)
        </div>
        <div className="sub">
          {Object.entries(d.files).map(([c, n]) => `${NAMES[c]}: ${n}`).join(" · ")} files · healthy range from the {d.settings.k_neighbours} normal files with the
          closest shaft speed · settings chosen on the <b>tune</b> half, results below on the <b>test</b> half · exported {d.generated}
        </div>
      </div>

      <div className="grid2">
        <div className="panel">
          <div className="panel-title bold">Detection on the TEST half</div>
          <table className="tbl wide">
            <thead>
              <tr><th>Class</th><th>Files</th><th>Envelope: fault / normal correct</th><th>Envelope: right fault type</th><th>Entropy: correct</th></tr>
            </thead>
            <tbody>
              {CLASSES.filter((c) => d.test[c]).map((c) => {
                const v = d.test[c];
                return (
                  <tr key={c}>
                    <td style={{ color: COLORS[c], fontWeight: 700 }}>{NAMES[c]}</td>
                    <td>{v.n}</td>
                    <td><Bar v={v.env_correct} /></td>
                    <td>{c === "normal" ? "–" : <Bar v={v.env_type_correct} />}</td>
                    <td><Bar v={v.ent_correct} /></td>
                  </tr>
                );
              })}
            </tbody>
          </table>
          <p className="muted small">For normal files "correct" means no false alarm. Rule: a fault is declared when its z-score against the healthy range is ≥ {d.settings.z_clear} (entropy: |z| ≥ {d.settings.z_clear}).</p>
        </div>

        <div className="panel">
          <div className="panel-title bold">Settings chosen automatically (tune half)</div>
          <table className="tbl wide">
            <thead><tr><th>Detector</th><th>Sensor</th><th>Band / measure</th><th>Tune hit</th><th>Tune false alarm</th></tr></thead>
            <tbody>
              {Object.entries(E).map(([c, e]) => (
                <tr key={c}>
                  <td style={{ color: COLORS[c], fontWeight: 700 }}>{NAMES[c]}</td>
                  <td>{e.sensor}</td>
                  <td>{e.band[0]}–{e.band[1]} Hz, {e.line} lines</td>
                  <td>{pct(e.tune_hit)}</td>
                  <td>{pct(e.tune_false_alarm)}</td>
                </tr>
              ))}
              <tr>
                <td style={{ fontWeight: 700 }}>Entropy (any fault)</td>
                <td>{N.sensor}</td>
                <td>{N.name}</td>
                <td>{pct(N.tune_hit)}</td>
                <td>{pct(N.tune_false_alarm)}</td>
              </tr>
            </tbody>
          </table>

          <div className="panel-title bold" style={{ marginTop: 12 }}>Confusion matrix — envelope method, test half</div>
          <table className="tbl wide">
            <thead><tr><th>true ↓ / found →</th>{CLASSES.map((c) => <th key={c}>{NAMES[c]}</th>)}</tr></thead>
            <tbody>
              {CLASSES.map((a) => (
                <tr key={a}>
                  <td style={{ color: COLORS[a], fontWeight: 700 }}>{NAMES[a]}</td>
                  {CLASSES.map((b) => (
                    <td key={b} style={{ fontWeight: a === b ? 700 : 400, background: a === b ? "#eaf6ea" : d.confusion_test[a][b] ? "#fdeaea" : undefined }}>
                      {d.confusion_test[a][b]}
                    </td>
                  ))}
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>

      <div className="panel">
        <div className="panel-title bold">
          z-score of the{" "}
          <select value={zKey} onChange={(e) => setZKey(e.target.value)}>
            {Object.keys(E).map((c) => <option key={c} value={c}>{NAMES[c]} detector</option>)}
          </select>{" "}
          against shaft speed (test files) — above the red line = fault declared
        </div>
        <ResponsiveContainer width="100%" height={320}>
          <ScatterChart margin={{ top: 10, right: 20, bottom: 24, left: 8 }}>
            <CartesianGrid strokeOpacity={0.3} />
            <XAxis dataKey="fr" type="number" name="shaft speed" unit=" Hz" domain={["dataMin - 1", "dataMax + 1"]} label={{ value: "Shaft speed (Hz)", position: "insideBottom", offset: -12 }} />
            <YAxis dataKey="z" type="number" name="z" width={50} />
            <ZAxis range={[36, 36]} />
            <Tooltip formatter={(v) => fmt(v, 2)} cursor={{ strokeDasharray: "3 3" }} />
            <ReferenceLine y={d.settings.z_clear} stroke={RED} strokeDasharray="5 3" />
            <Legend verticalAlign="top" height={24} />
            {CLASSES.map((c) => <Scatter key={c} name={NAMES[c]} data={scatter[c]} fill={COLORS[c]} fillOpacity={0.75} />)}
          </ScatterChart>
        </ResponsiveContainer>
      </div>

      <div className="panel">
        <div className="panel-title bold">By bearing and imbalance weight (test half)</div>
        <table className="tbl wide">
          <thead><tr><th>Bearing</th><th>Weight</th><th>Fault</th><th>Files</th><th>Envelope detected</th><th>Right type</th><th>Entropy detected</th></tr></thead>
          <tbody>
            {d.by_group_test.map((g, i) => (
              <tr key={i}>
                <td>{g.bearing}</td><td>{g.weight}</td>
                <td style={{ color: COLORS[g.cls], fontWeight: 600 }}>{NAMES[g.cls]}</td>
                <td>{g.n}</td><td>{pct(g.env_detected)}</td><td>{pct(g.env_type_correct)}</td><td>{pct(g.ent_detected)}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
      <p className="howto">
        Fault lines from the MaFaulDa bearing factors (× shaft speed): FTF {d.settings.fault_factors.FTF}, BSF {d.settings.fault_factors.BSF},
        BPFO {d.settings.fault_factors.BPFO}, BPFI {d.settings.fault_factors.BPFI}. Shaft speed from the tachometer of every file.
        Analysis functions from <code>{d.script}</code>.
      </p>
    </div>
  );
}

function Bar({ v }) {
  const c = v >= 0.9 ? GREEN : v >= 0.6 ? "#c98a00" : RED;
  return (
    <div style={{ display: "flex", alignItems: "center", gap: 8 }}>
      <div style={{ width: 120, height: 10, background: "#eee", borderRadius: 5, overflow: "hidden" }}>
        <div style={{ width: `${Math.round(v * 100)}%`, height: "100%", background: c }} />
      </div>
      <b style={{ color: c }}>{pct(v)}</b>
    </div>
  );
}
