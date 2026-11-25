import React, { useEffect, useMemo, useRef, useState } from "react";
import Window1 from "./Window1.jsx";
import Window2 from "./Window2.jsx";
import Window3 from "./Window3.jsx";
import Overview from "./Overview.jsx";
import MaFaulda from "./MaFaulda.jsx";
import { clock, GREEN, RED } from "./common.js";

const TABS = [
  ["w1", "1 · Live monitor"],
  ["w2", "2 · Entropy"],
  ["w3", "3 · VMD"],
  ["overview", "Overview"],
  ["mafaulda", "MaFaulDa dataset"],
];
const DATA = "./data/";

export default function App() {
  const [index, setIndex] = useState(null);
  const [error, setError] = useState(null);
  const [slug, setSlug] = useState(null);
  const [summary, setSummary] = useState(null);
  const [winNo, setWinNo] = useState(1);
  const [win, setWin] = useState(null);
  const [tab, setTab] = useState("w1");
  const [playing, setPlaying] = useState(false);
  const [speed, setSpeed] = useState(2); // seconds per 10 s window
  const cache = useRef({});

  // list of recordings
  useEffect(() => {
    fetch(DATA + "index.json")
      .then((r) => (r.ok ? r.json() : Promise.reject(new Error(`index.json: HTTP ${r.status}`))))
      .then((ix) => {
        setIndex(ix);
        const pick = ix.recordings.find((r) => r.label === "Faulty") || ix.recordings[0];
        if (pick) setSlug(pick.slug);
      })
      .catch((e) => setError(String(e)));
  }, []);

  // summary of the chosen recording
  useEffect(() => {
    if (!slug) return;
    setSummary(null);
    setWin(null);
    setPlaying(false);
    fetch(`${DATA}${slug}/summary.json`)
      .then((r) => r.json())
      .then((s) => {
        setSummary(s);
        setWinNo(s.windows.length);
      })
      .catch((e) => setError(String(e)));
  }, [slug]);

  // one 10 s window (cached)
  useEffect(() => {
    if (!summary) return;
    const key = `${slug}/${winNo}`;
    if (cache.current[key]) {
      setWin(cache.current[key]);
      return;
    }
    let alive = true;
    fetch(`${DATA}${slug}/w${String(winNo).padStart(2, "0")}.json`)
      .then((r) => r.json())
      .then((w) => {
        cache.current[key] = w;
        if (alive) setWin(w);
      })
      .catch((e) => setError(String(e)));
    return () => {
      alive = false;
    };
  }, [summary, slug, winNo]);

  // playback: step through the windows like the Python replay
  useEffect(() => {
    if (!playing || !summary) return;
    const id = setInterval(() => {
      setWinNo((n) => {
        if (n >= summary.windows.length) {
          setPlaying(false);
          return n;
        }
        return n + 1;
      });
    }, speed * 1000);
    return () => clearInterval(id);
  }, [playing, speed, summary]);

  // keyboard: SPACE play/pause, arrows step, R restart
  useEffect(() => {
    const onKey = (e) => {
      if (!summary || e.target.tagName === "SELECT" || e.target.tagName === "INPUT") return;
      if (e.code === "Space") {
        e.preventDefault();
        setPlaying((p) => !p);
      } else if (e.key === "ArrowRight") setWinNo((n) => Math.min(summary.windows.length, n + 1));
      else if (e.key === "ArrowLeft") setWinNo((n) => Math.max(1, n - 1));
      else if (e.key === "r" || e.key === "R") {
        setWinNo(1);
        setPlaying(true);
      }
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [summary]);

  const nWin = summary ? summary.windows.length : 0;
  const strip = useMemo(() => (summary ? summary.windows : []), [summary]);

  if (error)
    return (
      <div className="page">
        <div className="card error">
          <h2>Could not load the dashboard data</h2>
          <p>{error}</p>
          <p>
            Run <code>python export_dashboard_data.py</code> in the <code>dahsboard</code> folder first. It writes the
            results to <code>dashboard-app/public/data/</code>.
          </p>
        </div>
      </div>
    );
  if (!index) return <div className="page loading">Loading…</div>;

  return (
    <div className="page">
      <header className="topbar">
        <div className="brand">
          <div className="title">Bearing Fault Dashboard</div>
          <div className="subtitle">
            Fault feature extraction in vibration signals using entropy measures · results from{" "}
            <code>{index.script}</code> · exported {index.generated}
          </div>
        </div>
        <div className="controls">
          <label>
            Recording
            <select value={slug || ""} onChange={(e) => setSlug(e.target.value)}>
              {index.recordings.map((r) => (
                <option key={r.slug} value={r.slug}>
                  {r.slug} ({r.label})
                </option>
              ))}
            </select>
          </label>
          <button className="btn" onClick={() => setPlaying((p) => !p)} disabled={!summary}>
            {playing ? "❚❚ Pause" : "▶ Play"}
          </button>
          <button
            className="btn ghost"
            onClick={() => {
              setWinNo(1);
              setPlaying(true);
            }}
            disabled={!summary}
          >
            ↺ Restart
          </button>
          <label>
            Speed
            <select value={speed} onChange={(e) => setSpeed(Number(e.target.value))}>
              <option value={10}>x1 (real time)</option>
              <option value={2}>x5</option>
              <option value={1}>x10</option>
              <option value={0.5}>x20</option>
            </select>
          </label>
        </div>
      </header>

      {summary && (
        <section className="timeline card">
          <div className="timeline-head">
            <span>
              <b>Window {winNo}</b> / {nWin} · {clock((winNo - 1) * summary.windowS)}–{clock(winNo * summary.windowS)} of{" "}
              {clock(summary.durationS)} · Group {win ? win.group : "…"}/{summary.groups}
            </span>
            <span className="hint">SPACE play/pause · ← → step · R restart · click a square to jump</span>
          </div>
          <input
            type="range"
            min={1}
            max={nWin}
            value={winNo}
            onChange={(e) => {
              setPlaying(false);
              setWinNo(Number(e.target.value));
            }}
          />
          <div className="strip">
            {strip.map((w) => (
              <button
                key={w.index}
                title={`#${w.index}: ${w.verdict} · entropy ${w.entVerdict}`}
                className={"sq" + (w.index === winNo ? " on" : "")}
                onClick={() => {
                  setPlaying(false);
                  setWinNo(w.index);
                }}
              >
                <span style={{ background: w.verdict === "HEALTHY" ? GREEN : RED }} />
                <span style={{ background: w.entVerdict === "HEALTHY" ? GREEN : RED }} />
              </button>
            ))}
          </div>
          <div className="legend-small">top half = spectrum verdict (window 1) · bottom half = entropy verdict (window 2)</div>
        </section>
      )}

      <nav className="tabs">
        {TABS.map(([k, label]) => (
          <button key={k} className={tab === k ? "tab on" : "tab"} onClick={() => setTab(k)}>
            {label}
          </button>
        ))}
      </nav>

      <main>
        {tab === "mafaulda" ? (
          <MaFaulda />
        ) : tab === "overview" ? (
          <Overview index={index} onOpen={(s) => { setSlug(s); setTab("w1"); }} />
        ) : !summary || !win || win.index !== winNo ? (
          <div className="card loading">Loading window {winNo}…</div>
        ) : tab === "w1" ? (
          <Window1 s={summary} w={win} />
        ) : tab === "w2" ? (
          <Window2 s={summary} w={win} />
        ) : (
          <Window3 s={summary} w={win} />
        )}
      </main>
      <footer className="foot">
        Ammar Na'im bin Borhan · UiTM · FYP · Teensy 4.1 + ADXL326, Fs = {summary ? summary.fs : 2000} Hz
      </footer>
    </div>
  );
}
