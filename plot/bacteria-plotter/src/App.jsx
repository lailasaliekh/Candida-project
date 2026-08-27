import { useState, useRef, useCallback } from "react";

// ── Parser ────────────────────────────────────────────────────────────────────
function parseDat(text) {
  const lines = text.trim().split(/\r?\n/);
  if (lines.length < 2) throw new Error("File too short — need a header + data rows.");

  const headers = lines[0].trim().split(/\s+/);
  const idx = (name) => {
    const i = headers.indexOf(name);
    if (i === -1) throw new Error(`Missing column: "${name}"`);
    return i;
  };

  const cols = {
    cell_type:  idx("cell_type"),
    cell_id:    idx("cell_id"),
    length:     idx("length"),
    radius:     idx("radius"),
    pos_x:      idx("pos_x"),
    pos_y:      idx("pos_y"),
    pos_z:      idx("pos_z"),
    ori_x:      idx("ori_x"),
    ori_y:      idx("ori_y"),
    ori_z:      idx("ori_z"),
    lower_link: idx("lower_link"),
    upper_link: idx("upper_link"),
  };

  const parseLink = (v) => (v === "None" || v === "" || v === "-") ? null : parseInt(v, 10);

  return lines.slice(1).map((line, i) => {
    const f = line.trim().split(/\s+/);
    if (f.length < headers.length) throw new Error(`Row ${i + 2}: expected ${headers.length} fields, got ${f.length}.`);
    return {
      cell_type:  f[cols.cell_type],
      cell_id:    parseInt(f[cols.cell_id], 10),
      length:     parseFloat(f[cols.length]),
      radius:     parseFloat(f[cols.radius]),
      pos_x:      parseFloat(f[cols.pos_x]),
      pos_y:      parseFloat(f[cols.pos_y]),
      pos_z:      parseFloat(f[cols.pos_z]),
      ori_x:      parseFloat(f[cols.ori_x]),
      ori_y:      parseFloat(f[cols.ori_y]),
      ori_z:      parseFloat(f[cols.ori_z]),
      lower_link: parseLink(f[cols.lower_link]),
      upper_link: parseLink(f[cols.upper_link]),
    };
  });
}

// ── Plot ──────────────────────────────────────────────────────────────────────
function BacteriaPlot({ cells, filename }) {
  const W = 700, H = 520, SCALE = 30;

  const cellMap = Object.fromEntries(cells.map(c => [c.cell_id, c]));

  const chains = (() => {
    const pairs = [], seen = new Set();
    cells.forEach(c => {
      [c.upper_link, c.lower_link].forEach(link => {
        if (link && cellMap[link]) {
          const key = [c.cell_id, link].sort().join("-");
          if (!seen.has(key)) { seen.add(key); pairs.push([c.cell_id, link]); }
        }
      });
    });
    return pairs;
  })();

  const xs = cells.map(c => c.pos_x);
  const ys = cells.map(c => c.pos_y);
  const xMin = Math.min(...xs), xMax = Math.max(...xs);
  const yMin = Math.min(...ys), yMax = Math.max(...ys);
  const cx = (xMin + xMax) / 2;
  const cy = (yMin + yMax) / 2;

  // Auto-fit: zoom so entire colony (+ padding) fills the viewport
  const PAD = 2;
  const dataW = (xMax - xMin) + PAD * 2 + 2;
  const dataH = (yMax - yMin) + PAD * 2 + 2;
  const fitZoom = Math.min((W / dataW) / SCALE, (H / dataH) / SCALE, 4);

  const [hovered, setHovered]   = useState(null);
  const [pan, setPan]           = useState({ x: 0, y: 0 });
  const [zoom, setZoom]         = useState(fitZoom);
  const [dragging, setDragging] = useState(false);
  const lastMouse = useRef(null);

  const toSVG = (wx, wy) => ({
    x: W / 2 + (wx - cx) * SCALE * zoom + pan.x,
    y: H / 2 - (wy - cy) * SCALE * zoom + pan.y,
  });

  const chainType = (cell) =>
    [cell.lower_link, cell.upper_link].filter(Boolean).length === 0 ? "unchained" : "chained";

  const TYPE_STYLE = {
    unchained: { body:"url(#gradUnchained)", bodyHL:"url(#gradUnchainedHL)", stroke:"#ff9f43", glow:"#ff9f43", label:"#ffe0b2" },
    chained:   { body:"url(#gradChained)",   bodyHL:"url(#gradChainedHL)",   stroke:"#0af5a0", glow:"#0af5a0", label:"#c0ffe8" },
  };

  const capsule = (cell, rPx, highlighted) => {
    const { pos_x, pos_y, ori_x, ori_y, length } = cell;
    const { x, y } = toSVG(pos_x, pos_y);
    const cellLen  = length * SCALE * zoom;
    const angle    = Math.atan2(-ori_y, ori_x) * (180 / Math.PI);
    const ts       = TYPE_STYLE[chainType(cell)];
    return (
      <g key={cell.cell_id}
         transform={`translate(${x},${y}) rotate(${angle})`}
         style={{ cursor:"pointer" }}
         onMouseEnter={() => setHovered(cell.cell_id)}
         onMouseLeave={() => setHovered(null)}>
        <rect x={-cellLen/2} y={-rPx-3} width={cellLen} height={(rPx+3)*2}
              rx={rPx+3} ry={rPx+3} fill="none"
              stroke={ts.glow} strokeWidth={highlighted ? 4 : 2}
              opacity={highlighted ? 0.9 : 0.25} filter="url(#glow)" />
        <rect x={-cellLen/2} y={-rPx} width={cellLen} height={rPx*2}
              rx={rPx} ry={rPx}
              fill={highlighted ? ts.bodyHL : ts.body}
              stroke={ts.stroke} strokeWidth={highlighted ? 1.5 : 0.8} opacity={0.95} />
        {(highlighted || zoom > 1.5) && (
          <text x={0} y={-rPx-5} textAnchor="middle" fontSize={11/zoom}
                fill={ts.label} fontFamily="'Courier New',monospace"
                transform={`rotate(${-angle})`}
                style={{ pointerEvents:"none", userSelect:"none" }}>
            #{cell.cell_id}
          </text>
        )}
      </g>
    );
  };

  const handleWheel = (e) => {
    e.preventDefault();
    setZoom(z => Math.min(8, Math.max(0.2, z * (e.deltaY < 0 ? 1.1 : 0.9))));
  };
  const handleMouseDown = (e) => { setDragging(true); lastMouse.current = { x:e.clientX, y:e.clientY }; };
  const handleMouseMove = (e) => {
    if (!dragging) return;
    setPan(p => ({ x: p.x + e.clientX - lastMouse.current.x, y: p.y + e.clientY - lastMouse.current.y }));
    lastMouse.current = { x:e.clientX, y:e.clientY };
  };
  const handleMouseUp = () => setDragging(false);

  const rPx = 0.5 * SCALE * zoom;
  const hoveredCell = hovered ? cellMap[hovered] : null;

  // Dynamic grid
  const gridStep = 2;
  const gxMin = Math.floor((cx - W/(2*SCALE*zoom))/gridStep)*gridStep;
  const gxMax = Math.ceil( (cx + W/(2*SCALE*zoom))/gridStep)*gridStep;
  const gyMin = Math.floor((cy - H/(2*SCALE*zoom))/gridStep)*gridStep;
  const gyMax = Math.ceil( (cy + H/(2*SCALE*zoom))/gridStep)*gridStep;
  const gridXs = [], gridYs = [];
  for (let v = gxMin; v <= gxMax; v += gridStep) gridXs.push(v);
  for (let v = gyMin; v <= gyMax; v += gridStep) gridYs.push(v);

  const unchainedCount = cells.filter(c => chainType(c) === "unchained").length;
  const chainedCount   = cells.length - unchainedCount;

  return (
    <div style={{ display:"flex", flexDirection:"column", alignItems:"center", gap:"14px", width:"100%" }}>
      <div style={{ textAlign:"center" }}>
        <div style={{ color:"#0af5a0", fontSize:"10px", letterSpacing:"0.3em", textTransform:"uppercase", marginBottom:"3px" }}>
          {filename} · {cells.length} cells
        </div>
        <h1 style={{ color:"#c0ffe8", fontSize:"20px", fontWeight:400, margin:0, letterSpacing:"0.1em" }}>
          Bacterial Colony — XY Plane
        </h1>
      </div>

      <div style={{ position:"relative" }}>
        <svg width={W} height={H}
             style={{ border:"1px solid #0af5a030", borderRadius:"8px",
                      background:"radial-gradient(ellipse at center,#061a10 0%,#020a06 100%)",
                      cursor: dragging ? "grabbing" : "grab", display:"block" }}
             onWheel={handleWheel} onMouseDown={handleMouseDown}
             onMouseMove={handleMouseMove} onMouseUp={handleMouseUp} onMouseLeave={handleMouseUp}>
          <defs>
            <filter id="glow" x="-50%" y="-50%" width="200%" height="200%">
              <feGaussianBlur stdDeviation="3" result="coloredBlur"/>
              <feMerge><feMergeNode in="coloredBlur"/><feMergeNode in="SourceGraphic"/></feMerge>
            </filter>
            <linearGradient id="gradUnchained" x1="0" y1="0" x2="0" y2="1">
              <stop offset="0%"   stopColor="#ffcc70" stopOpacity="0.95"/>
              <stop offset="50%"  stopColor="#e08020" stopOpacity="0.88"/>
              <stop offset="100%" stopColor="#7a3800" stopOpacity="0.9"/>
            </linearGradient>
            <linearGradient id="gradUnchainedHL" x1="0" y1="0" x2="0" y2="1">
              <stop offset="0%"   stopColor="#ffe8a0" stopOpacity="1"/>
              <stop offset="50%"  stopColor="#ffaa44" stopOpacity="0.95"/>
              <stop offset="100%" stopColor="#a04800" stopOpacity="0.9"/>
            </linearGradient>
            <linearGradient id="gradChained" x1="0" y1="0" x2="0" y2="1">
              <stop offset="0%"   stopColor="#2affb0" stopOpacity="0.95"/>
              <stop offset="50%"  stopColor="#0dcc80" stopOpacity="0.85"/>
              <stop offset="100%" stopColor="#085c38" stopOpacity="0.9"/>
            </linearGradient>
            <linearGradient id="gradChainedHL" x1="0" y1="0" x2="0" y2="1">
              <stop offset="0%"   stopColor="#80ffdc" stopOpacity="1"/>
              <stop offset="50%"  stopColor="#00ffa0" stopOpacity="0.95"/>
              <stop offset="100%" stopColor="#00804d" stopOpacity="0.9"/>
            </linearGradient>
          </defs>

          {/* Grid */}
          {gridXs.map(gx => { const { x } = toSVG(gx, cy); return <line key={`gx${gx}`} x1={x} y1={0} x2={x} y2={H} stroke="#0af5a010" strokeWidth={1}/>; })}
          {gridYs.map(gy => { const { y } = toSVG(cx, gy); return <line key={`gy${gy}`} x1={0} y1={y} x2={W} y2={y} stroke="#0af5a010" strokeWidth={1}/>; })}
          {gridXs.filter((_,i) => i%2===0).map(gx => { const {x} = toSVG(gx,cy); return <text key={`lx${gx}`} x={x} y={H-4} textAnchor="middle" fontSize={8} fill="#0af5a040" fontFamily="'Courier New',monospace">{gx.toFixed(0)}</text>; })}
          {gridYs.filter((_,i) => i%2===0).map(gy => { const {y} = toSVG(cx,gy); return <text key={`ly${gy}`} x={4} y={y+3} fontSize={8} fill="#0af5a040" fontFamily="'Courier New',monospace">{gy.toFixed(0)}</text>; })}

          {/* Chain links */}
          {chains.map(([idA, idB]) => {
            const a = cellMap[idA], b = cellMap[idB];
            if (!a || !b) return null;
            const sa = toSVG(a.pos_x, a.pos_y), sb = toSVG(b.pos_x, b.pos_y);
            const isHov = hovered === idA || hovered === idB;
            return <line key={`${idA}-${idB}`} x1={sa.x} y1={sa.y} x2={sb.x} y2={sb.y}
                         stroke={isHov ? "#ffffff" : "#54a0ff55"}
                         strokeWidth={isHov ? 2 : 1} strokeDasharray={isHov ? "none" : "3,3"}/>;
          })}

          {/* Cells */}
          {cells.map(cell => capsule(cell, rPx, hovered === cell.cell_id))}

          <text x={W-8} y={H/2+14} textAnchor="end" fontSize={10} fill="#0af5a060" fontFamily="'Courier New',monospace">x (µm)</text>
          <text x={12}  y={14}      fontSize={10}    fill="#0af5a060" fontFamily="'Courier New',monospace">y (µm)</text>
        </svg>

        {hoveredCell && (
          <div style={{ position:"absolute", bottom:"12px", left:"12px",
                        background:"#071a10ee", border:"1px solid #0af5a050",
                        borderRadius:"6px", padding:"10px 14px", color:"#c0ffe8",
                        fontSize:"11px", lineHeight:"1.8", backdropFilter:"blur(4px)",
                        pointerEvents:"none" }}>
            <div style={{ color:"#0af5a0", marginBottom:"4px", fontWeight:"bold" }}>Cell #{hoveredCell.cell_id}</div>
            <div>Type: <span style={{color:"#fff"}}>{hoveredCell.cell_type}</span></div>
            <div>Length: <span style={{color:"#fff"}}>{hoveredCell.length.toFixed(3)} µm</span></div>
            <div>Radius: <span style={{color:"#fff"}}>{hoveredCell.radius} µm</span></div>
            <div>Position: <span style={{color:"#fff"}}>({hoveredCell.pos_x.toFixed(3)}, {hoveredCell.pos_y.toFixed(3)}, {hoveredCell.pos_z.toFixed(3)})</span></div>
            <div>Orientation: <span style={{color:"#fff"}}>({hoveredCell.ori_x.toFixed(3)}, {hoveredCell.ori_y.toFixed(3)})</span></div>
            <div>Links: <span style={{color:"#fff"}}>↓{hoveredCell.lower_link ?? "—"} ↑{hoveredCell.upper_link ?? "—"}</span></div>
          </div>
        )}
      </div>

      <div style={{ display:"flex", gap:"20px", alignItems:"center", fontSize:"10px",
                    letterSpacing:"0.1em", flexWrap:"wrap", justifyContent:"center" }}>
        {[
          { label:"Unchained", color:"#ff9f43", count: unchainedCount },
          { label:"Chained",   color:"#0af5a0", count: chainedCount   },
        ].map(({ label, color, count }) => (
          <div key={label} style={{ display:"flex", alignItems:"center", gap:"7px", color: color+"cc" }}>
            <svg width={28} height={12}>
              <rect x={0} y={2} width={28} height={8} rx={4} fill={color} fillOpacity={0.6} stroke={color} strokeWidth={0.8}/>
            </svg>
            {label} <span style={{ color:"#ffffff60" }}>×{count}</span>
          </div>
        ))}
        <button onClick={() => { setZoom(fitZoom); setPan({ x:0, y:0 }); }}
                style={{ background:"none", border:"1px solid #0af5a030", color:"#0af5a060",
                         borderRadius:"5px", padding:"3px 10px", cursor:"pointer",
                         fontSize:"10px", letterSpacing:"0.08em" }}>
          ⊡ Fit view
        </button>
        <span style={{ color:"#ffffff30" }}>· Scroll to zoom · Drag to pan · Hover for details</span>
      </div>
    </div>
  );
}

// ── Drop zone ─────────────────────────────────────────────────────────────────
function DropZone({ onLoad }) {
  const [dragOver, setDragOver] = useState(false);
  const [error, setError]       = useState(null);
  const inputRef = useRef(null);

  const readFile = useCallback((file) => {
    if (!file) return;
    setError(null);
    const reader = new FileReader();
    reader.onload = (e) => {
      try {
        const cells = parseDat(e.target.result);
        onLoad(cells, file.name);
      } catch (err) {
        setError(err.message);
      }
    };
    reader.readAsText(file);
  }, [onLoad]);

  const onDrop = (e) => {
    e.preventDefault();
    setDragOver(false);
    readFile(e.dataTransfer.files[0]);
  };

  return (
    <div style={{ display:"flex", flexDirection:"column", alignItems:"center", justifyContent:"center",
                  minHeight:"60vh", gap:"24px" }}>
      <div style={{ textAlign:"center" }}>
        <div style={{ color:"#0af5a0", fontSize:"11px", letterSpacing:"0.35em", textTransform:"uppercase", marginBottom:"6px" }}>
          Biofilm Visualiser
        </div>
        <h1 style={{ color:"#c0ffe8", fontSize:"26px", fontWeight:400, margin:0, letterSpacing:"0.1em" }}>
          Bacterial Colony Plotter
        </h1>
      </div>

      <div
        onClick={() => inputRef.current.click()}
        onDragOver={(e) => { e.preventDefault(); setDragOver(true); }}
        onDragLeave={() => setDragOver(false)}
        onDrop={onDrop}
        style={{
          width:"380px", height:"180px",
          border:`2px dashed ${dragOver ? "#0af5a0" : "#0af5a040"}`,
          borderRadius:"12px",
          background: dragOver ? "#0af5a012" : "#061a1060",
          display:"flex", flexDirection:"column", alignItems:"center", justifyContent:"center",
          gap:"12px", cursor:"pointer", transition:"all 0.2s ease",
        }}>
        <svg width={40} height={40} viewBox="0 0 24 24" fill="none"
             stroke={dragOver ? "#0af5a0" : "#0af5a060"} strokeWidth={1.5}>
          <path d="M4 16v2a2 2 0 002 2h12a2 2 0 002-2v-2M12 4v12M8 8l4-4 4 4" strokeLinecap="round" strokeLinejoin="round"/>
        </svg>
        <div style={{ color: dragOver ? "#0af5a0" : "#0af5a080", fontSize:"13px", textAlign:"center", lineHeight:"1.6" }}>
          Drop a <span style={{ color:"#c0ffe8" }}>biofilm*.dat</span> file here<br/>
          <span style={{ fontSize:"11px" }}>or click to browse</span>
        </div>
        <input ref={inputRef} type="file" accept=".dat,.txt"
               style={{ display:"none" }}
               onChange={(e) => readFile(e.target.files[0])} />
      </div>

      {error && (
        <div style={{ color:"#ff6b6b", fontSize:"12px", background:"#2a080830",
                      border:"1px solid #ff6b6b40", borderRadius:"6px",
                      padding:"10px 16px", maxWidth:"380px", lineHeight:"1.6" }}>
          ⚠ {error}
        </div>
      )}

      <div style={{ color:"#0af5a025", fontSize:"10px", maxWidth:"400px", textAlign:"center", lineHeight:"2" }}>
        Expected tab-separated columns:<br/>
        cell_type · cell_id · length · radius · pos_x · pos_y · pos_z<br/>
        ori_x · ori_y · ori_z · lower_link · upper_link
      </div>
    </div>
  );
}

// ── Root ──────────────────────────────────────────────────────────────────────
export default function App() {
  const [cells, setCells]       = useState(null);
  const [filename, setFilename] = useState("");

  return (
    <div style={{
      background:"#050d0a", minHeight:"100vh",
      display:"flex", flexDirection:"column",
      alignItems:"center", justifyContent:"center",
      fontFamily:"'Courier New', monospace", padding:"24px", gap:"16px",
    }}>
      {cells ? (
        <>
          <BacteriaPlot cells={cells} filename={filename} />
          <button onClick={() => setCells(null)}
                  style={{ background:"none", border:"1px solid #0af5a030", color:"#0af5a060",
                           borderRadius:"6px", padding:"6px 16px", cursor:"pointer",
                           fontSize:"11px", letterSpacing:"0.1em", marginTop:"4px" }}>
            ← Load another file
          </button>
        </>
      ) : (
        <DropZone onLoad={(c, n) => { setCells(c); setFilename(n); }} />
      )}
    </div>
  );
}
