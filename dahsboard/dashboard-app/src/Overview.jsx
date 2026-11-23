import React from "react";
import { GREEN, RED } from "./common.js";

export default function Overview({ index, onOpen }) {
  return (
    <div className="window">
      <div className="panel">
        <div className="panel-title bold">All recordings — how many 10 s windows each method calls FAULTY</div>
        <table className="tbl wide">
          <thead>
            <tr>
              <th>Recording</th><th>Label</th><th>Windows</th>
              <th>Spectrum verdict (window 1)</th><th>Entropy verdict (window 2)</th><th></th>
            </tr>
          </thead>
          <tbody>
            {index.recordings.map((r) => {
              const expectFault = r.label === "Faulty";
              const ok1 = expectFault ? r.faulty === r.windows : r.faulty === 0;
              const ok2 = expectFault ? r.entFaulty === r.windows : r.entFaulty === 0;
              return (
                <tr key={r.slug}>
                  <td><b>{r.slug}</b><div className="muted small">{r.name}.xlsx</div></td>
                  <td>{r.label}</td>
                  <td>{r.windows}</td>
                  <td style={{ color: ok1 ? GREEN : RED, fontWeight: 600 }}>{r.faulty} / {r.windows} FAULTY</td>
                  <td style={{ color: ok2 ? GREEN : RED, fontWeight: 600 }}>{r.entFaulty} / {r.windows} FAULTY</td>
                  <td><button className="btn ghost" onClick={() => onOpen(r.slug)}>Open</button></td>
                </tr>
              );
            })}
          </tbody>
        </table>
        <p className="muted small">
          Green = the method agrees with the recording's label in every window, red = at least one window disagrees.
          Results come from <code>{index.script}</code>, exported {index.generated}. To refresh after changing the
          Python code, run <code>python export_dashboard_data.py</code> again.
        </p>
      </div>
    </div>
  );
}
