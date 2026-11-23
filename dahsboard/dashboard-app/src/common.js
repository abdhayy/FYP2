// Shared colours and helpers. Colours follow the matplotlib "tab" palette used by the Python script.
export const FAULTS = ["BPFO", "BPFI", "BSF", "FTF"];
export const FAULT_COLORS = { BPFO: "#d62728", BPFI: "#ff7f0e", BSF: "#9467bd", FTF: "#2ca02c", SHAFT: "#7f7f7f" };
export const RED = "#d62728";
export const GREEN = "#008c00";
export const ORANGE = "#ff8c00";

export const fmt = (v, d = 1) => (v === null || v === undefined || Number.isNaN(v) ? "–" : Number(v).toFixed(d));
export const signed = (v, d = 1) => (v === null || v === undefined ? "–" : (v >= 0 ? "+" : "") + Number(v).toFixed(d));
export const clock = (s) => `${Math.floor(s / 60)}:${String(Math.round(s % 60)).padStart(2, "0")}`;

// zip parallel arrays into Recharts rows: zipRows({x: [...], y: [...]}) -> [{x, y}, ...]
export function zipRows(cols) {
  const keys = Object.keys(cols);
  const n = Math.min(...keys.map((k) => cols[k].length));
  const out = new Array(n);
  for (let i = 0; i < n; i++) {
    const row = {};
    for (const k of keys) row[k] = cols[k][i];
    out[i] = row;
  }
  return out;
}

// lines of the bearing faults drawn on the plain FFT (same as the Python script: BPFO, BPFI, 2xBSF, FTF)
export function fftMarkers(b, known) {
  return [
    { key: "BPFO", f: b.BPFO, name: "BPFO" + (known.BPFO ? "" : "*") },
    { key: "BPFI", f: b.BPFI, name: "BPFI" + (known.BPFI ? "" : "*") },
    { key: "BSF", f: 2 * b.BSF, name: "2xBSF" },
    { key: "FTF", f: b.FTF, name: "FTF" },
  ];
}

// lines drawn on the envelope FFT (same list as mark_fault_lines in the Python script)
export function envMarkers(b) {
  const fr = b.fr;
  return [
    ["BPFO", b.BPFO, "BPFO", "solid"], ["BPFO", 2 * b.BPFO, "2xBPFO", "dash"], ["BPFO", 3 * b.BPFO, "3xBPFO", "dash"],
    ["BPFI", b.BPFI, "BPFI", "solid"], ["BPFI", b.BPFI - fr, "BPFI-fr", "dot"], ["BPFI", b.BPFI + fr, "BPFI+fr", "dot"],
    ["BSF", 2 * b.BSF, "2xBSF", "solid"], ["BSF", b.BSF, "BSF", "dash"],
    ["FTF", b.FTF, "FTF", "solid"], ["FTF", 2 * b.FTF, "2xFTF", "dash"],
  ].filter(([, f]) => f < 450).map(([fault, f, name, ls]) => ({ fault, f, name, ls }));
}

export const dashOf = (ls) => (ls === "dash" ? "6 3" : ls === "dot" ? "2 3" : undefined);
