const OPTION_TYPE = { CALL: "call", PUT: "put" };
const OPTION_STYLE = {
  EUROPEAN: "european",
  AMERICAN: "american",
  ASIAN: "asian",
  BARRIER: "barrier",
  LOOKBACK: "lookback",
  DIGITAL: "digital",
};
const BARRIER_TYPE = {
  UP_AND_OUT: "up_and_out",
  UP_AND_IN: "up_and_in",
  DOWN_AND_OUT: "down_and_out",
  DOWN_AND_IN: "down_and_in",
};

const DEFAULT_STATE = {
  underlying: "SPX",
  S: 100,
  K: 100,
  T: 1,
  r: 0.05,
  sigma: 0.2,
  q: 0,
  optionType: OPTION_TYPE.CALL,
  style: OPTION_STYLE.EUROPEAN,
  useHeston: false,
  barrier: 115,
  barrierType: BARRIER_TYPE.UP_AND_OUT,
  rebate: 0,
  avgType: "arithmetic",
  averagingFreq: 252,
  lookbackType: "floating",
  cashPayout: 10,
  sims: 20000,
  treeSteps: 250,
  outlook: "neutral",
  volView: "neutral",
  riskTolerance: "moderate",
};

const state = { ...DEFAULT_STATE };
const $ = (id) => document.getElementById(id);

function erf(x) {
  const sign = x < 0 ? -1 : 1;
  const a1 = 0.254829592;
  const a2 = -0.284496736;
  const a3 = 1.421413741;
  const a4 = -1.453152027;
  const a5 = 1.061405429;
  const p = 0.3275911;
  const t = 1 / (1 + p * Math.abs(x));
  const y = 1 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
  return sign * y;
}

function normCdf(x) {
  return 0.5 * (1 + erf(x / Math.sqrt(2)));
}

function normPdf(x) {
  return Math.exp(-0.5 * x * x) / Math.sqrt(2 * Math.PI);
}

function average(values) {
  return values.reduce((sum, value) => sum + value, 0) / values.length;
}

function stdev(values) {
  const mean = average(values);
  return Math.sqrt(average(values.map((value) => (value - mean) ** 2)));
}

function formatCurrency(value) {
  return `$${value.toFixed(2)}`;
}

function formatRange([low, high]) {
  return `${formatCurrency(low)} - ${formatCurrency(high)}`;
}

function randomNormal(seed) {
  let value = seed >>> 0;
  return () => {
    value = (1664525 * value + 1013904223) >>> 0;
    const u1 = Math.max(value / 4294967296, 1e-12);
    value = (1664525 * value + 1013904223) >>> 0;
    const u2 = value / 4294967296;
    return Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
  };
}

function createContract(overrides = {}) {
  const contract = {
    S: state.S,
    K: state.K,
    T: state.T,
    r: state.r,
    sigma: state.sigma,
    q: state.q,
    optionType: state.optionType,
    style: state.style,
    barrier: state.barrier,
    barrierType: state.barrierType,
    rebate: state.rebate,
    avgType: state.avgType,
    averagingFreq: state.averagingFreq,
    lookbackType: state.lookbackType,
    cashPayout: state.cashPayout,
    sims: state.sims,
    treeSteps: state.treeSteps,
    ...overrides,
  };
  if (contract.S <= 0 || contract.K <= 0 || contract.T <= 0 || contract.sigma <= 0) {
    throw new Error("Spot, strike, expiry, and volatility must be positive.");
  }
  return contract;
}

function isCall(contract) {
  return contract.optionType === OPTION_TYPE.CALL;
}

function moneyness(contract) {
  const ratio = contract.S / contract.K;
  if (ratio > 1.05) return isCall(contract) ? "ITM" : "OTM";
  if (ratio < 0.95) return isCall(contract) ? "OTM" : "ITM";
  return "ATM";
}

function d1d2(contract) {
  const d1 = (Math.log(contract.S / contract.K) +
    (contract.r - contract.q + 0.5 * contract.sigma ** 2) * contract.T) /
    (contract.sigma * Math.sqrt(contract.T));
  return { d1, d2: d1 - contract.sigma * Math.sqrt(contract.T) };
}

class BlackScholesModel {
  price(contract) {
    const { d1, d2 } = d1d2(contract);
    const discQ = Math.exp(-contract.q * contract.T);
    const discR = Math.exp(-contract.r * contract.T);
    return isCall(contract)
      ? contract.S * discQ * normCdf(d1) - contract.K * discR * normCdf(d2)
      : contract.K * discR * normCdf(-d2) - contract.S * discQ * normCdf(-d1);
  }

  greeks(contract) {
    const { d1, d2 } = d1d2(contract);
    const discQ = Math.exp(-contract.q * contract.T);
    const discR = Math.exp(-contract.r * contract.T);
    const sqrtT = Math.sqrt(contract.T);
    const pdf = normPdf(d1);
    const thetaCore = -(contract.S * discQ * pdf * contract.sigma) / (2 * sqrtT);
    return {
      delta: isCall(contract) ? discQ * normCdf(d1) : -discQ * normCdf(-d1),
      gamma: discQ * pdf / (contract.S * contract.sigma * sqrtT),
      theta: isCall(contract)
        ? (thetaCore - contract.r * contract.K * discR * normCdf(d2) + contract.q * contract.S * discQ * normCdf(d1)) / 365
        : (thetaCore + contract.r * contract.K * discR * normCdf(-d2) - contract.q * contract.S * discQ * normCdf(-d1)) / 365,
      vega: contract.S * discQ * pdf * sqrtT / 100,
      rho: isCall(contract)
        ? contract.K * contract.T * discR * normCdf(d2) / 100
        : -contract.K * contract.T * discR * normCdf(-d2) / 100,
      vanna: -discQ * pdf * d2 / contract.sigma,
      volga: contract.S * discQ * pdf * sqrtT * d1 * d2 / contract.sigma,
    };
  }
}

class BinomialTreeModel {
  price(contract) {
    const N = Math.max(25, Math.floor(contract.treeSteps));
    const dt = contract.T / N;
    const u = Math.exp(contract.sigma * Math.sqrt(dt));
    const d = 1 / u;
    const disc = Math.exp(-contract.r * dt);
    const p = (Math.exp((contract.r - contract.q) * dt) - d) / (u - d);
    const values = [];
    for (let j = 0; j <= N; j += 1) {
      const ST = contract.S * (u ** (N - j)) * (d ** j);
      values[j] = isCall(contract) ? Math.max(ST - contract.K, 0) : Math.max(contract.K - ST, 0);
    }
    for (let i = N - 1; i >= 0; i -= 1) {
      for (let j = 0; j <= i; j += 1) {
        const hold = disc * (p * values[j] + (1 - p) * values[j + 1]);
        if (contract.style === OPTION_STYLE.AMERICAN) {
          const SNode = contract.S * (u ** (i - j)) * (d ** j);
          const intrinsic = isCall(contract) ? Math.max(SNode - contract.K, 0) : Math.max(contract.K - SNode, 0);
          values[j] = Math.max(hold, intrinsic);
        } else {
          values[j] = hold;
        }
      }
    }
    return values[0];
  }
}

function payoffForPath(contract, path) {
  const ST = path[path.length - 1];
  if (contract.style === OPTION_STYLE.EUROPEAN || contract.style === OPTION_STYLE.AMERICAN) {
    return isCall(contract) ? Math.max(ST - contract.K, 0) : Math.max(contract.K - ST, 0);
  }
  if (contract.style === OPTION_STYLE.ASIAN) {
    const sample = path.slice(1);
    const avg = contract.avgType === "geometric"
      ? Math.exp(sample.reduce((sum, price) => sum + Math.log(price), 0) / sample.length)
      : average(sample);
    return isCall(contract) ? Math.max(avg - contract.K, 0) : Math.max(contract.K - avg, 0);
  }
  if (contract.style === OPTION_STYLE.BARRIER) {
    const raw = isCall(contract) ? Math.max(ST - contract.K, 0) : Math.max(contract.K - ST, 0);
    const top = Math.max(...path);
    const bottom = Math.min(...path);
    if (contract.barrierType === BARRIER_TYPE.UP_AND_OUT) return top >= contract.barrier ? contract.rebate : raw;
    if (contract.barrierType === BARRIER_TYPE.UP_AND_IN) return top >= contract.barrier ? raw : contract.rebate;
    if (contract.barrierType === BARRIER_TYPE.DOWN_AND_OUT) return bottom <= contract.barrier ? contract.rebate : raw;
    return bottom <= contract.barrier ? raw : contract.rebate;
  }
  if (contract.style === OPTION_STYLE.LOOKBACK) {
    if (contract.lookbackType === "floating") {
      return isCall(contract) ? ST - Math.min(...path) : Math.max(...path) - ST;
    }
    return isCall(contract)
      ? Math.max(Math.max(...path) - contract.K, 0)
      : Math.max(contract.K - Math.min(...path), 0);
  }
  return 0;
}

class MonteCarloModel {
  price(contract) {
    const sims = Math.max(4000, Math.floor(contract.sims));
    const steps = Math.max(25, Math.floor(contract.averagingFreq));
    const rng = randomNormal(42);
    const dt = contract.T / steps;
    const drift = (contract.r - contract.q - 0.5 * contract.sigma ** 2) * dt;
    const vol = contract.sigma * Math.sqrt(dt);
    const discounted = [];
    for (let i = 0; i < sims; i += 1) {
      let spot = contract.S;
      const path = [spot];
      for (let j = 0; j < steps; j += 1) {
        spot *= Math.exp(drift + vol * rng());
        path.push(spot);
      }
      discounted.push(Math.exp(-contract.r * contract.T) * payoffForPath(contract, path));
    }
    const price = average(discounted);
    const stderr = stdev(discounted) / Math.sqrt(sims);
    return { price, stderr, ci95: [price - 1.96 * stderr, price + 1.96 * stderr], sims };
  }
}

class HestonOverlay {
  price(contract, bsmPrice) {
    const term = Math.sqrt(contract.T) * 0.07;
    const skew = (contract.K / contract.S - 1) * (isCall(contract) ? -0.4 : 0.25);
    return Math.max(0, bsmPrice * (1 + 0.015 + term + skew));
  }
}

class GreeksEngine {
  numerical(contract, pricer) {
    const dS = contract.S * 0.01;
    const dSig = contract.sigma * 0.01;
    const dr = 0.0001;
    const oneDay = 1 / 365;
    const make = (overrides = {}) => createContract({ ...contract, ...overrides });
    const p0 = pricer(contract);
    const pUp = pricer(make({ S: contract.S + dS }));
    const pDn = pricer(make({ S: contract.S - dS }));
    const pSigUp = pricer(make({ sigma: contract.sigma + dSig }));
    const pSigDn = pricer(make({ sigma: Math.max(0.0001, contract.sigma - dSig) }));
    const pRate = pricer(make({ r: contract.r + dr }));
    const pTime = pricer(make({ T: Math.max(0.0001, contract.T - oneDay) }));
    return {
      delta: (pUp - pDn) / (2 * dS),
      gamma: (pUp - 2 * p0 + pDn) / (dS ** 2),
      vega: ((pSigUp - pSigDn) / (2 * dSig)) * 0.01,
      theta: (pTime - p0) / (-oneDay),
      rho: ((pRate - p0) / dr) * 0.01,
    };
  }

  payoffAtExpiry(contracts, multipliers, premiums) {
    const strikes = contracts.map((contract) => contract.K);
    const minSpot = Math.min(...strikes) * 0.6;
    const maxSpot = Math.max(...strikes) * 1.4;
    return Array.from({ length: 120 }, (_, index) => {
      const spot = minSpot + (index / 119) * (maxSpot - minSpot);
      const total = contracts.reduce((sum, contract, i) => {
        const payoff = contract.optionType === OPTION_TYPE.CALL
          ? Math.max(spot - contract.K, 0)
          : Math.max(contract.K - spot, 0);
        return sum + multipliers[i] * (payoff - premiums[i]);
      }, 0);
      return { x: spot, y: total };
    });
  }
}

class StrategyRecommender {
  constructor(bsm) {
    this.bsm = bsm;
    this.catalog = {
      long_call: { legs: [["call", 1, 0]], outlook: "Bullish", vol: "Long" },
      long_put: { legs: [["put", 1, 0]], outlook: "Bearish", vol: "Long" },
      bull_call_spread: { legs: [["call", 1, 0], ["call", -1, 5]], outlook: "Moderately Bullish", vol: "Neutral" },
      bear_put_spread: { legs: [["put", 1, 0], ["put", -1, -5]], outlook: "Moderately Bearish", vol: "Neutral" },
      long_straddle: { legs: [["call", 1, 0], ["put", 1, 0]], outlook: "Neutral", vol: "Long" },
      short_straddle: { legs: [["call", -1, 0], ["put", -1, 0]], outlook: "Neutral", vol: "Short" },
      long_strangle: { legs: [["call", 1, 5], ["put", 1, -5]], outlook: "Neutral", vol: "Long" },
      short_strangle: { legs: [["call", -1, 5], ["put", -1, -5]], outlook: "Neutral", vol: "Short" },
      iron_condor: { legs: [["put", -1, -5], ["put", 1, -10], ["call", -1, 5], ["call", 1, 10]], outlook: "Neutral", vol: "Short" },
      risk_reversal: { legs: [["call", 1, 5], ["put", -1, -5]], outlook: "Bullish", vol: "Neutral" },
    };
  }

  recommend(contract, outlook, volView, riskTolerance) {
    const scored = Object.entries(this.catalog).map(([name, meta]) => {
      let score = 0;
      if (meta.outlook.toLowerCase().includes(outlook)) score += 3;
      if (meta.vol.toLowerCase().includes(volView)) score += 2;
      if (riskTolerance === "moderate" && meta.legs.length === 2) score += 1;
      if (riskTolerance === "low" && meta.legs.length >= 3) score += 1;
      if (riskTolerance === "high" && meta.legs.length === 1) score += 1;
      return { name, score, ...meta };
    }).sort((a, b) => b.score - a.score).slice(0, 4);

    return scored.map((item) => {
      const contracts = [];
      const multipliers = [];
      const premiums = [];
      const legs = item.legs.map(([type, mult, offset]) => {
        const legContract = createContract({ K: contract.K + offset, optionType: type, style: OPTION_STYLE.EUROPEAN });
        const premium = this.bsm.price(legContract);
        contracts.push(legContract);
        multipliers.push(mult);
        premiums.push(premium);
        const side = mult > 0 ? "Long" : "Short";
        return `${side} ${type.toUpperCase()} ${legContract.K.toFixed(0)} @ ${premium.toFixed(2)}`;
      });
      const net = premiums.reduce((sum, premium, index) => sum + premium * multipliers[index], 0);
      return {
        strategy: item.name.replaceAll("_", " "),
        score: item.score,
        outlook: item.outlook,
        volView: item.vol,
        legs,
        contracts,
        multipliers,
        premiums,
        net,
      };
    });
  }
}

class UnifiedPricer {
  constructor() {
    this.bsm = new BlackScholesModel();
    this.binomial = new BinomialTreeModel();
    this.mc = new MonteCarloModel();
    this.heston = new HestonOverlay();
    this.greeks = new GreeksEngine();
  }

  digitalPrice(contract) {
    const { d2 } = d1d2(contract);
    const disc = Math.exp(-contract.r * contract.T);
    return disc * contract.cashPayout * (isCall(contract) ? normCdf(d2) : normCdf(-d2));
  }

  price(contract, useHeston = false) {
    const result = { modelUsed: "", price: 0, greeks: {}, notes: [] };
    if (contract.style === OPTION_STYLE.DIGITAL) {
      result.modelUsed = "BSM Digital";
      result.price = this.digitalPrice(contract);
      result.greeks = this.greeks.numerical(contract, (c) => this.digitalPrice(c));
      result.notes.push("Cash-or-nothing payout.");
      return result;
    }

    if (contract.style === OPTION_STYLE.EUROPEAN) {
      const bsmPrice = this.bsm.price(contract);
      result.modelUsed = "Black-Scholes";
      result.price = bsmPrice;
      result.greeks = this.bsm.greeks(contract);
      result.notes.push("Closed-form routing for vanilla European options.");
      if (useHeston) {
        result.hestonPrice = this.heston.price(contract, bsmPrice);
        result.smilePremium = result.hestonPrice - bsmPrice;
        result.modelUsed = "Black-Scholes + Heston Overlay";
        result.price = result.hestonPrice;
        result.notes.push(`Smile premium ${formatCurrency(result.smilePremium)}.`);
      }
      return result;
    }

    if (contract.style === OPTION_STYLE.AMERICAN) {
      const price = this.binomial.price(contract);
      const euro = this.bsm.price(createContract({ ...contract, style: OPTION_STYLE.EUROPEAN }));
      result.modelUsed = "Binomial Tree";
      result.price = price;
      result.earlyExercisePremium = price - euro;
      result.greeks = this.greeks.numerical(contract, (c) => this.binomial.price(c));
      result.notes.push(`Early exercise premium ${formatCurrency(result.earlyExercisePremium)}.`);
      return result;
    }

    if ([OPTION_STYLE.ASIAN, OPTION_STYLE.BARRIER, OPTION_STYLE.LOOKBACK].includes(contract.style)) {
      const mc = this.mc.price(contract);
      result.modelUsed = "Monte Carlo";
      result.price = mc.price;
      result.stderr = mc.stderr;
      result.ci95 = mc.ci95;
      result.greeks = this.greeks.numerical(contract, (c) => this.mc.price(c).price);
      result.notes.push(`95% confidence band ${formatRange(mc.ci95)}.`);
      return result;
    }

    throw new Error(`Unsupported style ${contract.style}`);
  }
}

function generateHeatmap(contract) {
  const maturities = [0.25, 0.5, 1, 1.5, 2];
  const strikes = [0.8, 0.9, 1, 1.1, 1.2].map((ratio) => contract.S * ratio);
  return maturities.flatMap((T) => strikes.map((K) => ({
    T,
    K,
    iv: (contract.sigma + Math.abs(K / contract.S - 1) * 0.18 + 0.02 * Math.sqrt(T)) * 100,
  })));
}

function generateSmile(contract) {
  return [-20, -10, -5, 0, 5, 10, 20].map((offset) => ({
    x: Math.max(1, contract.K + offset),
    y: (Math.max(0.05, contract.sigma + Math.abs(offset) / 200) * 100),
  }));
}

function addLog(message) {
  const entry = document.createElement("div");
  entry.className = "terminal-entry";
  const time = document.createElement("div");
  time.className = "time";
  time.textContent = new Intl.DateTimeFormat("en-US", { hour: "2-digit", minute: "2-digit", second: "2-digit" }).format(new Date());
  const body = document.createElement("div");
  body.textContent = message;
  entry.append(time, body);
  $("terminalLog").prepend(entry);
}

function renderDiagnostics(result, contract) {
  const items = [
    { label: "Model", value: result.modelUsed },
    { label: "Underlying", value: `${state.underlying} ${formatCurrency(contract.S)}` },
    { label: "Moneyness", value: moneyness(contract) },
    { label: "Style", value: contract.style },
  ];
  if (result.hestonPrice !== undefined) items.push({ label: "Heston Price", value: formatCurrency(result.hestonPrice) });
  if (result.earlyExercisePremium !== undefined) items.push({ label: "Early Exercise", value: formatCurrency(result.earlyExercisePremium) });
  if (result.ci95) items.push({ label: "95% Band", value: formatRange(result.ci95) });
  result.notes.forEach((note, index) => items.push({ label: `Note ${index + 1}`, value: note }));
  $("diagnostics").innerHTML = "";
  items.forEach((item) => {
    const row = document.createElement("div");
    row.className = "diagnostic-item";
    row.innerHTML = `<div class="diagnostic-label">${item.label}</div><div class="diagnostic-value">${item.value}</div>`;
    $("diagnostics").appendChild(row);
  });
}

function renderHeatmap(data) {
  const root = $("ivHeatmap");
  root.innerHTML = "";
  const values = data.map((item) => item.iv);
  const min = Math.min(...values);
  const max = Math.max(...values);
  data.forEach((item) => {
    const ratio = (item.iv - min) / Math.max(max - min, 1e-9);
    const hue = 200 - ratio * 160;
    const cell = document.createElement("div");
    cell.className = "heatmap-cell";
    cell.style.background = `linear-gradient(180deg, hsl(${hue} 90% 78%), hsl(${hue - 20} 95% 60%))`;
    cell.innerHTML = `<span>K ${item.K.toFixed(0)} | T ${item.T.toFixed(2)}</span><strong>${item.iv.toFixed(1)}%</strong>`;
    root.appendChild(cell);
  });
}

function renderLineChart(svg, points, title, xLabel, stroke) {
  const width = 640;
  const height = 280;
  const pad = 28;
  const xs = points.map((point) => point.x);
  const ys = points.map((point) => point.y);
  const xMin = Math.min(...xs);
  const xMax = Math.max(...xs);
  const yMin = Math.min(...ys);
  const yMax = Math.max(...ys);
  const sx = (value) => pad + ((value - xMin) / Math.max(xMax - xMin, 1e-9)) * (width - 2 * pad);
  const sy = (value) => height - pad - ((value - yMin) / Math.max(yMax - yMin, 1e-9)) * (height - 2 * pad);
  const zeroY = sy(Math.min(Math.max(0, yMin), yMax));
  const path = points.map((point, index) => `${index === 0 ? "M" : "L"} ${sx(point.x)} ${sy(point.y)}`).join(" ");
  svg.innerHTML = `
    <line x1="${pad}" y1="${zeroY}" x2="${width - pad}" y2="${zeroY}" stroke="rgba(255,255,255,0.18)" stroke-dasharray="5 5"></line>
    <line x1="${pad}" y1="${pad}" x2="${pad}" y2="${height - pad}" stroke="rgba(255,255,255,0.18)"></line>
    <path d="${path}" fill="none" stroke="${stroke}" stroke-width="4" stroke-linecap="round"></path>
    <text x="${pad}" y="18" fill="#80a0ba" font-family="Consolas" font-size="11">${title}</text>
    <text x="${width - pad - 100}" y="${height - 10}" fill="#80a0ba" font-family="Consolas" font-size="11">${xLabel}</text>
  `;
}

function syncInputs() {
  $("underlyingInput").value = state.underlying;
  $("spotInput").value = state.S;
  $("strikeInput").value = state.K;
  $("expiryInput").value = state.T;
  $("rateInput").value = state.r;
  $("volInput").value = state.sigma;
  $("dividendInput").value = state.q;
  $("optionTypeInput").value = state.optionType;
  $("styleInput").value = state.style;
  $("hestonToggleInput").value = state.useHeston ? "on" : "off";
  $("barrierInput").value = state.barrier;
  $("barrierTypeInput").value = state.barrierType;
  $("rebateInput").value = state.rebate;
  $("avgTypeInput").value = state.avgType;
  $("avgFreqInput").value = state.averagingFreq;
  $("lookbackTypeInput").value = state.lookbackType;
  $("cashPayoutInput").value = state.cashPayout;
  $("simsInput").value = state.sims;
  $("treeStepsInput").value = state.treeSteps;
  $("outlookInput").value = state.outlook;
  $("volViewInput").value = state.volView;
  $("riskInput").value = state.riskTolerance;
}

function collectFormState() {
  state.underlying = $("underlyingInput").value.trim() || "SPX";
  state.S = Number($("spotInput").value);
  state.K = Number($("strikeInput").value);
  state.T = Number($("expiryInput").value);
  state.r = Number($("rateInput").value);
  state.sigma = Number($("volInput").value);
  state.q = Number($("dividendInput").value);
  state.optionType = $("optionTypeInput").value;
  state.style = $("styleInput").value;
  state.useHeston = $("hestonToggleInput").value === "on";
  state.barrier = Number($("barrierInput").value);
  state.barrierType = $("barrierTypeInput").value;
  state.rebate = Number($("rebateInput").value);
  state.avgType = $("avgTypeInput").value;
  state.averagingFreq = Number($("avgFreqInput").value);
  state.lookbackType = $("lookbackTypeInput").value;
  state.cashPayout = Number($("cashPayoutInput").value);
  state.sims = Number($("simsInput").value);
  state.treeSteps = Number($("treeStepsInput").value);
  state.outlook = $("outlookInput").value;
  state.volView = $("volViewInput").value;
  state.riskTolerance = $("riskInput").value;
}

function renderStrategies(contract) {
  const strategies = new StrategyRecommender(new BlackScholesModel())
    .recommend(contract, state.outlook, state.volView, state.riskTolerance);
  $("strategyList").innerHTML = "";
  strategies.forEach((strategy) => {
    const card = document.createElement("div");
    card.className = "strategy-card";
    card.innerHTML = `
      <h3>${strategy.strategy}</h3>
      <div class="meta">
        <span>Score ${strategy.score}</span>
        <span>${strategy.outlook}</span>
        <span>${strategy.volView} vol</span>
        <span>Net ${formatCurrency(strategy.net)}</span>
      </div>
      <div class="legs">${strategy.legs.join("<br>")}</div>
    `;
    card.addEventListener("click", () => {
      const payoff = new GreeksEngine().payoffAtExpiry(strategy.contracts, strategy.multipliers, strategy.premiums);
      renderLineChart($("payoffChart"), payoff, `Payoff | ${strategy.strategy}`, "Underlying at Expiry", "#f2a93b");
      addLog(`Loaded ${strategy.strategy} payoff profile.`);
    });
    $("strategyList").appendChild(card);
  });
}

function updateDashboard(result, contract) {
  const greeks = result.greeks || {};
  $("heroSpot").textContent = contract.S.toFixed(2);
  $("heroPrice").textContent = result.price.toFixed(2);
  $("heroModel").textContent = result.modelUsed;
  $("heroSignal").textContent = state.outlook.charAt(0).toUpperCase() + state.outlook.slice(1);
  $("theoValue").textContent = formatCurrency(result.price);
  $("deltaValue").textContent = (greeks.delta || 0).toFixed(4);
  $("gammaValue").textContent = (greeks.gamma || 0).toFixed(4);
  $("vegaValue").textContent = (greeks.vega || 0).toFixed(4);
  $("thetaValue").textContent = (greeks.theta || 0).toFixed(4);
  $("rhoValue").textContent = (greeks.rho || 0).toFixed(4);
  $("priceAnnotation").textContent = result.modelUsed;
  $("volRegime").textContent = state.sigma > 0.35 ? "High Vol" : state.sigma < 0.18 ? "Low Vol" : "Balanced";
  $("moneynessBadge").textContent = moneyness(contract);
  $("confidenceBand").textContent = result.ci95 ? formatRange(result.ci95) : "n/a";
  renderDiagnostics(result, contract);
  renderHeatmap(generateHeatmap(contract));
  renderLineChart($("smileChart"), generateSmile(contract), "Synthetic IV Smile", "Strike", "#38d4ff");
  const singleLeg = createContract({ style: OPTION_STYLE.EUROPEAN });
  const premium = new BlackScholesModel().price(singleLeg);
  const payoff = new GreeksEngine().payoffAtExpiry([singleLeg], [1], [premium]);
  renderLineChart($("payoffChart"), payoff, "Single-Leg Payoff", "Underlying at Expiry", "#f2a93b");
  renderStrategies(contract);
}

function priceCurrentTicket(logIt = true) {
  try {
    collectFormState();
    const contract = createContract();
    const result = new UnifiedPricer().price(contract, state.useHeston);
    updateDashboard(result, contract);
    if (logIt) {
      addLog(`Priced ${state.underlying} ${state.optionType.toUpperCase()} ${state.style.toUpperCase()} at ${formatCurrency(result.price)} via ${result.modelUsed}.`);
    }
  } catch (error) {
    addLog(`Error: ${error.message}`);
  }
}

function runCommandString() {
  const text = $("commandInput").value.trim();
  const upper = text.toUpperCase();
  if (upper.includes("CALL")) state.optionType = OPTION_TYPE.CALL;
  if (upper.includes("PUT")) state.optionType = OPTION_TYPE.PUT;
  if (upper.includes("EUROPEAN")) state.style = OPTION_STYLE.EUROPEAN;
  if (upper.includes("AMERICAN")) state.style = OPTION_STYLE.AMERICAN;
  if (upper.includes("ASIAN")) state.style = OPTION_STYLE.ASIAN;
  if (upper.includes("BARRIER")) state.style = OPTION_STYLE.BARRIER;
  if (upper.includes("LOOKBACK")) state.style = OPTION_STYLE.LOOKBACK;
  if (upper.includes("DIGITAL")) state.style = OPTION_STYLE.DIGITAL;
  state.useHeston = upper.includes("HESTON");
  syncInputs();
  priceCurrentTicket();
  addLog(`Executed terminal command: ${text}`);
}

function updateClock() {
  $("deskClock").textContent = new Intl.DateTimeFormat("en-US", {
    hour: "2-digit",
    minute: "2-digit",
    second: "2-digit",
  }).format(new Date());
}

$("pricingForm").addEventListener("submit", (event) => {
  event.preventDefault();
  priceCurrentTicket();
});

$("runCommandBtn").addEventListener("click", runCommandString);
$("commandInput").addEventListener("keydown", (event) => {
  if (event.key === "Enter") runCommandString();
});

$("resetDefaultsBtn").addEventListener("click", () => {
  Object.assign(state, DEFAULT_STATE);
  syncInputs();
  priceCurrentTicket(false);
  addLog("Reset pricing ticket to default desk values.");
});

["outlookInput", "volViewInput", "riskInput"].forEach((id) => {
  $(id).addEventListener("change", () => {
    collectFormState();
    renderStrategies(createContract());
    addLog("Updated strategy scanner filters.");
  });
});

setInterval(updateClock, 1000);
updateClock();
syncInputs();
priceCurrentTicket(false);
addLog("VOL<GO> terminal initialized.");
