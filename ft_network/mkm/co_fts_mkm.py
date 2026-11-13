#!/usr/bin/env python3
"""钴基费托合成的集总微观动力学模型 (lumped MKM)。

单一位点类型 (代表 Co 台阶位 / B5 位系综), 平均场近似, 每个吸附物种占 1 个位点。

网络 (全部基元步骤, 括号内为参数名):
  CO(g) + *      ⇌ CO*                 (E_CO, w_CO: CO 吸附能及其覆盖度相关的侧向排斥)
  H2(g) + 2*     ⇌ 2H*                 (E_H2, w_H: 每个 H2 的解离吸附能)
  CO* + *        ⇌ C* + O*             (dir)   CO 直接解离
  CO* + H*       ⇌ HCO* + *            (hco)   氢助活化第一步
  HCO* + *       ⇌ CH* + O*            (hcod)  HCO 的 C–O 断裂
  C* + H*        ⇌ CH* + *             (ch)
  CH* + H*       ⇌ CH2* + *            (ch2)
  CH2* + H*      ⇌ CH3* + *            (ch3)
  CH3* + H*      →  CH4(g) + 2*        (ch4)   甲烷化
  O* + H*        ⇌ OH* + *             (oh)
  OH* + H*       →  H2O(g) + 2*        (h2o)
  CO* + O*       →  CO2(g) + 2*        (co2)
  R_n* + CH2*    →  R_{n+1}* + *       (grow)  链增长 (R_1 = CH3*, R_n = CnH2n+1* 烷基)
  R_n* + H*      →  CnH2n+2(g) + 2*    (par)   链终止: 烷烃 (n ≥ 2)
  R_n* + *       ⇌  E_n* + H*          (beta)  β-H 消除 → 吸附 α-烯烃 (n ≥ 2)
  E_n*           →  CnH2n(g) + *       (des)   烯烃脱附, 脱附能随 n 线性增加: E_des0 + c_des·(n−2)

“烯烃再吸附”用 E_n* 的再加氢 (beta 逆反应) 与脱附的竞争来表示; 脱附能随链长增加 → 长链烯烃更容易被
再加氢并继续增长, 这对应 Iglesia / Kuipers 解释的 α 随链长增大与烯烃/烷烃比随链长下降。

速率常数: 表面步骤 k = (kB T/h)·exp(−Ea/kBT), 逆反应 Ea − ΔE (忽略吸附物熵变)。
吸附: 碰撞理论 (粘附系数 S0); 脱附: 由吸附平衡常数 (气相熵全部损失) 求出。
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

KB_EV = 8.617333e-5           # eV/K
KB = 1.380649e-23             # J/K
H_PL = 6.62607015e-34         # J s
AMU = 1.66053907e-27
A_SITE = 5.44e-20             # m², Co(0001) 每个表面原子的面积
P0 = 1e5                      # Pa


def s_gas(name, T):
    """气相标准熵 (eV/K), S298 + Cp·ln(T/298) 近似。"""
    s298, cp = {"CO": (197.7, 29.1), "H2": (130.7, 28.8)}[name]
    return (s298 + cp * math.log(T / 298.15)) / 96485.0


DEFAULT = dict(
    # 吸附 (eV, 正值 = 放热的结合能)
    E_CO=1.50, w_CO=1.30, S0_CO=1.0,
    E_H2=0.90, w_H=0.30, S0_H2=0.1,
    # 表面步骤: (Ea, ΔE) eV
    dir=(1.10, -0.40), hco=(1.00, 0.70), hcod=(0.90, -0.90),
    ch=(0.70, -0.20), ch2=(0.65, 0.10), ch3=(0.60, -0.15), ch4=(0.95, -0.10),
    oh=(1.20, 0.10), h2o=(1.00, -0.10), co2=(1.30, -0.20),
    grow=(0.80, -0.50), par=(1.00, -0.40), beta=(0.75, 0.25),
    E_des0=0.55, c_des=0.035,
)


@dataclass
class Result:
    T: float
    P: float
    H2CO: float
    theta: dict
    tof_co: float
    s_ch4: float
    s_co2: float
    paraffin: np.ndarray          # 每个 n 的生成速率 (s⁻¹), 下标 n
    olefin: np.ndarray
    resid: float
    extra: dict = field(default_factory=dict)

    def carbon_sel(self):
        """各碳数产物的碳基选择性 (按 CO 消耗计)。"""
        n = np.arange(len(self.paraffin))
        return n * (self.paraffin + self.olefin) / self.tof_co

    def weight_frac(self, nmin=1):
        n = np.arange(len(self.paraffin))
        w = n * (self.paraffin + self.olefin)
        w[1] = self.paraffin[1]                 # n = 1 只有 CH4
        w[:nmin] = 0
        return w / w.sum()

    def alpha(self, lo=10, hi=20):
        n = np.arange(len(self.paraffin))
        tot = self.paraffin + self.olefin
        sel = (n >= lo) & (n <= hi)
        return math.exp(np.polyfit(n[sel], np.log(tot[sel]), 1)[0])

    def op_ratio(self, n):
        return self.olefin[n] / self.paraffin[n]


class CoMKM:
    NFIX = ["CO", "H", "O", "OH", "C", "CH", "CH2", "HCO"]

    def __init__(self, params=None, N=40):
        self.p = dict(DEFAULT)
        if params:
            self.p.update(params)
        self.N = N
        self.idx = {s: i for i, s in enumerate(self.NFIX)}
        self.iR = len(self.NFIX)                 # R_1 … R_N
        self.iE = self.iR + N                    # E_2 … E_N
        self.nvar = self.iE + (N - 1)
        self._last = None

    # ---------------------------------------------------------------- 常数
    def _consts(self, T, P, h2co, pp=None):
        p, kT = self.p, KB_EV * T
        pre = KB * T / H_PL
        f = lambda Ea: pre * np.exp(-np.asarray(Ea, float) / kT)                     # noqa: E731
        k = {}
        for name in ["dir", "hco", "hcod", "ch", "ch2", "ch3", "ch4", "oh", "h2o", "co2", "grow", "par", "beta"]:
            Ea, dE = p[name]
            k[name] = (f(Ea), f(Ea - dE))
        if pp is None:
            PCO = P * 1e5 / (1 + h2co)
            PH2 = P * 1e5 * h2co / (1 + h2co)
        else:                                   # 直接给分压 (bar)
            PCO, PH2 = pp[0] * 1e5, pp[1] * 1e5
        flux = lambda m: A_SITE / math.sqrt(2 * math.pi * m * AMU * KB * T)   # noqa: E731  (每 Pa)
        ka_co, ka_h2 = p["S0_CO"] * flux(28.0), p["S0_H2"] * flux(2.016)
        k["co_ads"] = ka_co * PCO
        k["h2_ads"] = ka_h2 * PH2
        # 脱附常数 kd(θ) = ka·P0·exp(−(E_b(θ) − T·S_gas)/kT)
        k["co_des0"] = ka_co * P0 * math.exp(-(p["E_CO"] - T * s_gas("CO", T)) / kT)
        k["h2_des0"] = ka_h2 * P0 * math.exp(-(p["E_H2"] - T * s_gas("H2", T)) / kT)
        k["w_co"], k["w_h"], k["kT"] = p["w_CO"], p["w_H"], kT
        n = np.arange(2, self.N + 1)
        k["des"] = f(p["E_des0"] + p["c_des"] * (n - 2))
        self.k = k

    # ---------------------------------------------------------------- 速率
    def rates(self, th):
        k, i, N = self.k, self.idx, self.N
        CO, H, O, OH, C, CH, CH2, HCO = (th[i[s]] for s in self.NFIX)
        R = th[self.iR:self.iR + N]               # R[0] = R_1
        E = th[self.iE:self.iE + N - 1]           # E[0] = E_2
        free = 1.0 - th.sum()
        tot_cov = th.sum()
        kd_co = k["co_des0"] * math.exp(k["w_co"] * CO / k["kT"])
        kd_h2 = k["h2_des0"] * math.exp(2 * k["w_h"] * H / k["kT"])
        r = {}
        r["co_ads"] = k["co_ads"] * free - kd_co * CO
        r["h2_ads"] = k["h2_ads"] * free ** 2 - kd_h2 * H ** 2
        r["dir"] = k["dir"][0] * CO * free - k["dir"][1] * C * O
        r["hco"] = k["hco"][0] * CO * H - k["hco"][1] * HCO * free
        r["hcod"] = k["hcod"][0] * HCO * free - k["hcod"][1] * CH * O
        r["ch"] = k["ch"][0] * C * H - k["ch"][1] * CH * free
        r["ch2"] = k["ch2"][0] * CH * H - k["ch2"][1] * CH2 * free
        r["ch3"] = k["ch3"][0] * CH2 * H - k["ch3"][1] * R[0] * free
        r["ch4"] = k["ch4"][0] * R[0] * H
        r["oh"] = k["oh"][0] * O * H - k["oh"][1] * OH * free
        r["h2o"] = k["h2o"][0] * OH * H
        r["co2"] = k["co2"][0] * CO * O
        g = k["grow"][0] * R[:-1] * CH2           # R_n → R_{n+1}, n = 1…N−1
        par = k["par"][0] * R[1:] * H             # n = 2…N
        beta = k["beta"][0] * R[1:] * free - k["beta"][1] * E * H
        des = k["des"] * E
        return r, g, par, beta, des, free, tot_cov

    def deriv(self, th):
        r, g, par, beta, des, free, _ = self.rates(th)
        i, N = self.idx, self.N
        d = np.zeros_like(th)
        d[i["CO"]] = r["co_ads"] - r["dir"] - r["hco"] - r["co2"]
        d[i["H"]] = (2 * r["h2_ads"] - r["hco"] - r["ch"] - r["ch2"] - r["ch3"] - r["ch4"] - r["oh"] - r["h2o"]
                     - par.sum() + beta.sum())
        d[i["O"]] = r["dir"] + r["hcod"] - r["oh"] - r["co2"]
        d[i["OH"]] = r["oh"] - r["h2o"]
        d[i["C"]] = r["dir"] - r["ch"]
        d[i["CH"]] = r["hcod"] + r["ch"] - r["ch2"]
        d[i["CH2"]] = r["ch2"] - r["ch3"] - g.sum()
        d[i["HCO"]] = r["hco"] - r["hcod"]
        dR = np.zeros(N)
        dR[0] = r["ch3"] - r["ch4"] - g[0]
        dR[1:] += g
        dR[1:-1] -= g[1:]
        dR[1:] -= par + beta
        d[self.iR:self.iR + N] = dR
        d[self.iE:self.iE + N - 1] = beta - des
        return d

    # ---------------------------------------------------------------- 稳态
    def solve(self, T=493.15, P=20.0, h2co=2.0, warm=True, pp=None):
        self._consts(T, P, h2co, pp)
        x0 = self._last if (warm and self._last is not None) else None
        sol = None
        if x0 is not None:
            sol = root(lambda x: self.deriv(x), x0, method="hybr", options=dict(xtol=1e-12))
            if not (sol.success and np.all(sol.x > -1e-12) and sol.x.sum() < 1):
                sol = None
        if sol is None:
            th0 = np.zeros(self.nvar)
            ivp = solve_ivp(lambda t, x: self.deriv(x), (0, 1e9), th0, method="LSODA", rtol=1e-9, atol=1e-16)
            x = np.clip(ivp.y[:, -1], 0, None)
            sol = root(lambda y: self.deriv(y), x, method="hybr", options=dict(xtol=1e-12))
            if not (sol.success and np.all(sol.x > -1e-12)):
                sol.x = x
        th = np.clip(sol.x, 0, None)
        self._last = th
        r, g, par, beta, des, free, _ = self.rates(th)
        Np = self.N + 1
        paraffin, olefin = np.zeros(Np), np.zeros(Np)
        paraffin[1] = r["ch4"]
        paraffin[2:] = par
        olefin[2:] = des
        tof = r["co_ads"]
        theta = {s: th[self.idx[s]] for s in self.NFIX}
        theta["R(CnH2n+1)"] = th[self.iR:self.iR + self.N].sum()
        theta["CH3"] = th[self.iR]
        theta["E(olefin)"] = th[self.iE:].sum()
        theta["*"] = free
        resid = float(np.abs(self.deriv(th)).max())
        return Result(T, P, h2co, theta, tof, r["ch4"] / tof if tof > 0 else float("nan"),
                      r["co2"] / tof if tof > 0 else float("nan"), paraffin, olefin, resid,
                      extra=dict(rates=r, grow=g.sum()))


def kinetics(m, T=493.15, P=20.0, h2co=2.0, dT=5.0, f=1.2):
    """表观活化能 (kJ/mol) 和 H2、CO 反应级数 (以 CO 消耗 TOF 计)。"""
    PCO, PH2 = P / (1 + h2co), P * h2co / (1 + h2co)
    base = m.solve(T, pp=(PCO, PH2)).tof_co
    hi = m.solve(T + dT, pp=(PCO, PH2)).tof_co
    lo = m.solve(T - dT, pp=(PCO, PH2)).tof_co
    Ea = -8.314e-3 * (math.log(hi) - math.log(lo)) / (1 / (T + dT) - 1 / (T - dT))
    nH = math.log(m.solve(T, pp=(PCO, PH2 * f)).tof_co / base) / math.log(f)
    nCO = math.log(m.solve(T, pp=(PCO * f, PH2)).tof_co / base) / math.log(f)
    m.solve(T, pp=(PCO, PH2))
    return Ea, nH, nCO


if __name__ == "__main__":
    import time
    m = CoMKM()
    t0 = time.time()
    res = m.solve()
    print(f"{time.time() - t0:.2f} s, resid {res.resid:.2e}")
    print("TOF", res.tof_co, "S_CH4", res.s_ch4, "alpha", res.alpha(), "O/P C3", res.op_ratio(3))
    print({k: f"{v:.3g}" for k, v in res.theta.items()})
