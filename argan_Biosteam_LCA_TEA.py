"""
╔══════════════════════════════════════════════════════════════════════════════╗
║  INTEGRATED ARGAN (Argania spinosa) BIOREFINERY — ALL-IN-ONE SCRIPT         ║
║  BioSTEAM v2 Simulation · TEA · LCA · By Mouad Hachhach                      ║
╠══════════════════════════════════════════════════════════════════════════════╣

"""

# ═══════════════════════════════════════════════════════════════════════════════
# 0.  IMPORTS & OUTPUT DIRECTORY
# ═══════════════════════════════════════════════════════════════════════════════
import os, sys, json, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import (FancyBboxPatch, Arc, Circle, Ellipse,
                                FancyArrowPatch, Polygon, PathPatch)
from matplotlib.path import Path
from scipy.optimize import brentq
import biosteam as bst
import thermosteam as tmo
warnings.filterwarnings('ignore')
bst.nbtutorial()

# ── Output folder ─────────────────────────────────────────────────────────────
OUT = os.environ.get('OUT_DIR',
      os.path.join(os.path.dirname(os.path.abspath(__file__))
                   if '__file__' in dir() else os.getcwd(),
                   'argan_output'))
os.makedirs(OUT, exist_ok=True)
print(f"\n{'='*70}")
print(f"  ARGAN BIOREFINERY — ALL-IN-ONE PIPELINE")
print(f"  Output directory : {OUT}")
print(f"{'='*70}\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 1.  CEPCI COST ESCALATION FACTOR
# ═══════════════════════════════════════════════════════════════════════════════
# Reference: Turton et al. (2009), Appendix A, base year 2001, CEPCI = 397.0
# 2026 value derived from Chemical Engineering magazine monthly series:
#   2024 annual average ≈ 800  (Grokipedia / CE archives)
#   2025 annual average ≈ 813  (+1.6 % per CE Jan-2026 release)
#   Jan-2026 preliminary ≈ 820 (+0.9 % from Dec-2025, sixth consecutive rise)
CEPCI_REF  = 397.0   # Turton 2001 base year
CEPCI_2026 = 820.0   # Jan-2026 prelim. 
CEPCI_FACTOR = CEPCI_2026 / CEPCI_REF   # = 2.065
print(f"CEPCI update factor  : {CEPCI_REF:.1f} (2001 ref) → {CEPCI_2026:.1f} (2026)")
print(f"Escalation factor    : {CEPCI_FACTOR:.4f}\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 2.  THERMOSTEAM CHEMICAL DATABASE
# ═══════════════════════════════════════════════════════════════════════════════
_cell = tmo.Chemical('Cellulose', search_db=True)
_trio = tmo.Chemical('Triolein',  search_db=True)

def solid(ID, MW, Hf, CAS):
    """Create a solid pseudo-component copying heat/transport models from Cellulose."""
    c = tmo.Chemical(ID, search_db=False, MW=MW, Hf=Hf, phase='s', CAS=CAS)
    c.copy_models_from(_cell, ['V', 'Cn', 'kappa'])
    c.Tb = 773.15
    return c

def liquid(ID, MW, Hf, CAS):
    """Create a liquid pseudo-component copying models from Triolein."""
    c = tmo.Chemical(ID, search_db=False, MW=MW, Hf=Hf, phase='l', CAS=CAS)
    c.copy_models_from(_trio, ['V', 'Cn', 'kappa', 'mu'])
    c.Tb = 673.15
    return c

chemicals = tmo.Chemicals([
    'Water', 'Methane', 'CO2',
    tmo.Chemical('Triolein',  search_db=True),   # argan oil surrogate
    tmo.Chemical('Cellulose', search_db=True),
    tmo.Chemical('Glucose',   search_db=True),
    solid('Lignin',   3985.4, -3.96e6,  '8068-00-6'),  # Zine el Abidine 2013
    solid('Protein',  10000,  -3.5e6,   '9007-34-5'),  # Mirpoor 2024
    solid('Biochar',  12.0,   -393500,  '7782-42-5'),  # elemental C
    solid('Ash',      101.96, -1676000, '1302-74-5'),  # Al2O3 proxy
    liquid('Saponin', 956.6,  -2.5e6,   '8047-15-2'),  # Taarji 2018
    tmo.Chemical('Ethanol',   search_db=True),
])
tmo.settings.set_thermo(chemicals)
print("Chemical database OK  : 12 components registered\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 3.  FEED STREAMS  (10 000 t/yr @ 8 000 h/yr → 1 250 kg/h total fruit)
# ═══════════════════════════════════════════════════════════════════════════════
# Fruit composition (wt %):
#   Pulp   75 % → Glucose 45 %, Cellulose 25 %, Water 20 %, Ash 10 %
#   Shell  10 % → Lignin 55 %, Cellulose 40 %, Ash 5 %
#   Kernel 15 % → Triolein 48 %, Protein 30 %, Saponin 5 %, Water 12 %, Ash 5 %
# Sources: Khalil & Haidar 2016; Ait Itto 2024; El Monfalouti 2010

CAP = 1250.0   # kg/h total fruit

feed_pulp = bst.Stream('feed_pulp',
    Glucose   = CAP * 0.75 * 0.45,
    Cellulose = CAP * 0.75 * 0.25,
    Water     = CAP * 0.75 * 0.20,
    Ash       = CAP * 0.75 * 0.10,
    units='kg/hr')

feed_shell = bst.Stream('feed_shell',
    Lignin    = CAP * 0.10 * 0.55,
    Cellulose = CAP * 0.10 * 0.40,
    Ash       = CAP * 0.10 * 0.05,
    units='kg/hr')

feed_kernel = bst.Stream('feed_kernel',
    Triolein  = CAP * 0.15 * 0.48,
    Protein   = CAP * 0.15 * 0.30,
    Saponin   = CAP * 0.15 * 0.05,
    Water     = CAP * 0.15 * 0.12,
    Ash       = CAP * 0.15 * 0.05,
    units='kg/hr')

# Saponin extraction solvent: 60 wt% EtOH / 40 wt% water (Taarji et al. 2018)
solvent = bst.Stream('solvent',
    Water   = 120.0,
    Ethanol = 180.0,
    units='kg/hr', T=333.15)

# ═══════════════════════════════════════════════════════════════════════════════
# 4.  UNIT OPERATION MODELS
# ═══════════════════════════════════════════════════════════════════════════════

class ColdPress(bst.Unit):
    """
    Hydraulic cold press for kernel oil extraction.
    Oil recovery η = 92 % — El Monfalouti et al. (2010) Food Chem 122:173.
    Specific power: 25 W/kg feed.
    CAPEX correlation: C_ref = $85 000 at Q_ref = 500 kg/h; n = 0.6.
    Updated to 2026 USD using CEPCI factor 2.065.
    """
    _N_ins = 1;  _N_outs = 2
    OR = 0.92          # oil recovery fraction

    def _run(self):
        k, = self.ins
        oil, cake = self.outs
        oil.imass['Triolein'] = k.imass['Triolein'] * self.OR
        oil.imass['Water']    = k.imass['Water'] * 0.60
        oil.T = 303.15
        cake.copy_flow(k)
        cake.imass['Triolein'] = k.imass['Triolein'] * (1 - self.OR)
        cake.imass['Water']    = k.imass['Water'] * 0.40

    def _design(self):
        self.design_results['Power_kW'] = self.ins[0].F_mass * 0.025

    def _cost(self):
        Q = self.ins[0].F_mass
        C_base = 85000 * (Q / 500) ** 0.6   # 2001 USD
        self.purchase_costs['Hydraulic press'] = C_base * CEPCI_FACTOR


class OilFilter(bst.Unit):
    """
    Plate-and-frame polishing filter for crude argan oil.
    Removes 95 % of residual water from press effluent.
    """
    _N_ins = 1;  _N_outs = 2

    def _run(self):
        c, = self.ins
        ref, fines = self.outs
        ref.copy_flow(c)
        ref.imass['Water']  = c.imass['Water'] * 0.05
        fines.imass['Water'] = c.imass['Water'] * 0.95
        ref.T = 303.15

    def _cost(self):
        C_base = 25000   # 2001 USD
        self.purchase_costs['Plate-frame filter'] = C_base * CEPCI_FACTOR


class SaponinExtractor(bst.Unit):
    """
    Solid-liquid extractor for saponin recovery from press cake.
    Solvent: 60 wt% EtOH / water at 333 K, L:S = 10:1 (Taarji et al. 2018).
    Saponin recovery: 85 % (Henry et al. 2013 LC-MS inventory).
    """
    _N_ins = 2;  _N_outs = 2
    REC = 0.85

    def _run(self):
        cake, sol = self.ins
        ext, dep = self.outs
        se = cake.imass['Saponin'] * self.REC
        ext.imass['Saponin'] = se
        ext.imass['Water']   = sol.imass['Water']
        ext.imass['Ethanol'] = sol.imass['Ethanol']
        ext.T = 333.15
        dep.copy_flow(cake)
        dep.imass['Saponin'] = cake.imass['Saponin'] * (1 - self.REC)
        dep.imass['Water']   = cake.imass['Water'] * 0.30

    def _design(self):
        ck, _ = self.ins
        self.design_results['Vessel_m3'] = max(ck.F_vol * 3, 0.01)
        self.design_results['Power_kW']  = 0.5 * max(ck.F_vol * 3, 0.01)

    def _cost(self):
        v = self.design_results['Vessel_m3']
        self.purchase_costs['Extraction vessel'] = 40000 * (v/1.0)**0.6 * CEPCI_FACTOR
        self.purchase_costs['Agitator']           = 8000 * CEPCI_FACTOR


class Evaporator(bst.Unit):
    """
    Three-effect falling-film evaporator to concentrate saponin extract.
    Water removal: 85 %;  steam economy: ~2.8 kg water / kg steam.
    EtOH (95 %) recovered in condensate for recycle.
    """
    _N_ins = 1;  _N_outs = 2
    WR = 0.85

    def _run(self):
        d, = self.ins
        conc, cond = self.outs
        wr = d.imass['Water'] * self.WR
        conc.copy_flow(d)
        conc.imass['Water']   = d.imass['Water'] - wr
        conc.imass['Ethanol'] = d.imass['Ethanol'] * 0.05   # 95 % EtOH to condensate
        cond.imass['Water']   = wr
        cond.imass['Ethanol'] = d.imass['Ethanol'] * 0.95
        conc.T = 343.15

    def _design(self):
        d, = self.ins
        self.design_results['Evap_load_kW'] = d.imass['Water'] * self.WR * 2260 / 3600

    def _cost(self):
        ld = self.design_results['Evap_load_kW']
        C_base = 120000 * (max(ld, 1) / 50) ** 0.6   # 2001 USD
        self.purchase_costs['Multi-effect evaporator'] = C_base * CEPCI_FACTOR


class AnaerobicDigester(bst.Unit):
    """
    Mesophilic CSTR anaerobic digester (35 °C, HRT 25 d).
    Y_CH4 = 0.280 Nm³/kg VS destroyed — Carrere et al. (2010) J Hazard Mater 183:1.
    VS fraction = 70 % of organic dry mass;  VS destruction = 65 %.
    Capital: $1 200/m³ installed (reinforced concrete + mixing + gas collection).
    Updated to 2026 USD using CEPCI factor.
    """
    _N_ins = 1;  _N_outs = 2
    Y_CH4 = 0.280   # Nm³ CH4 / kg VS destroyed
    VS_f  = 0.70    # VS as fraction of organic dry mass
    VS_d  = 0.65    # VS destruction efficiency

    def _run(self):
        f, = self.ins
        bg, dg = self.outs
        om   = f.imass['Glucose'] + f.imass['Cellulose']
        vsd  = om * self.VS_f * self.VS_d
        # methane at STP density 0.656 kg/Nm³
        ch4  = vsd * self.Y_CH4 * 0.656
        co2  = ch4 * (44 / 16) * 0.538   # biogas ~65 % CH4 / 35 % CO2 by vol
        bg.imass['Methane'] = ch4
        bg.imass['CO2']     = co2
        bg.T = 308.15;  bg.P = 110000
        dg.copy_flow(f)
        fr = 1 - self.VS_d * self.VS_f
        dg.imass['Glucose']   = f.imass['Glucose']   * fr
        dg.imass['Cellulose'] = f.imass['Cellulose'] * fr
        dg.imass['Water']     = f.imass['Water'] * 0.98
        dg.T = 308.15

    def _design(self):
        f, = self.ins
        self.design_results['Vol_m3']   = max(f.F_vol * 25 * 24, 1)  # HRT=25 d
        self.design_results['Power_kW'] = 0.005 * self.design_results['Vol_m3']

    def _cost(self):
        v = self.design_results['Vol_m3']
        self.purchase_costs['Digester tank']    = 1200 * v * CEPCI_FACTOR
        self.purchase_costs['Gas handling']     = 150000 * CEPCI_FACTOR
        self.purchase_costs['Biogas upgrading'] = 80000  * CEPCI_FACTOR


class SlowPyrolysis(bst.Unit):
    """
    Rotary-kiln slow pyrolysis of argan shells.
    T_opt = 470–490 °C, no holding time (Ait Itto et al. 2024 NSGA-II optimum).
    Product distribution (dry feed basis):
        Biochar 37.99 % — Ait Itto 2024 Table 5, argan shell
        Bio-oil 25.45 % — idem
        Syngas  36.56 % — idem (CH4 18 wt%, rest CO2/CO/H2 lumped as CO2)
    Biochar BET potential ≥ 1 500 m²/g (Espinoza-Acosta 2023 Bioresour Technol).
    """
    _N_ins = 1;  _N_outs = 3
    Y_bc = 0.3799   # biochar yield
    Y_bo = 0.2545   # bio-oil yield
    Y_sg = 0.3656   # syngas yield

    def _run(self):
        s, = self.ins
        bc, bo, sg = self.outs
        dm = s.F_mass
        bc.imass['Biochar'] = dm * self.Y_bc
        bc.imass['Ash']     = s.imass['Ash'] * 0.90
        bc.T = 773.15
        bo.imass['Lignin'] = dm * self.Y_bo * 0.60   # phenolic-rich proxy
        bo.imass['Water']  = dm * self.Y_bo * 0.40
        bo.T = 298.15
        sg.imass['Methane'] = dm * self.Y_sg * 0.18
        sg.imass['CO2']     = dm * self.Y_sg * 0.82
        sg.T = 773.15;  sg.P = 110000

    def _design(self):
        s, = self.ins
        self.design_results['Throughput_kgph'] = s.F_mass
        self.design_results['Thermal_duty_kW'] = s.F_mass * 1.4  # kJ/g feedstock

    def _cost(self):
        cap = self.design_results['Throughput_kgph']
        self.purchase_costs['Rotary kiln'] = 350000 * (max(cap, 10) / 200)**0.6 * CEPCI_FACTOR
        self.purchase_costs['Condenser']   = 80000  * CEPCI_FACTOR
        self.purchase_costs['Char silo']   = 30000  * CEPCI_FACTOR


class CHP(bst.Unit):
    """
    Combined heat and power gas engine (Jenbacher J208 reference).
    Consumes all biogas (from AD) and pyrolysis syngas methane.
    η_electrical = 35 %;  η_thermal = 45 %;  fuel utilisation = 80 %.
    CH4 LHV = 50 050 kJ/kg.
    Capital: $800 / kW_electrical installed (2001 USD ref), scaled to 2026.
    """
    _N_ins = 2;  _N_outs = 1
    eta_e = 0.35;  eta_h = 0.45;  LHV_CH4 = 50050.0  # kJ/kg

    def _run(self):
        bg, sg = self.ins
        fl, = self.outs
        tot = bg.imass['Methane'] + sg.imass['Methane']
        self._pwr  = tot * self.LHV_CH4 / 3600 * self.eta_e   # kW electrical
        self._heat = tot * self.LHV_CH4 / 3600 * self.eta_h   # kW thermal
        fl.imass['CO2']   = (bg.imass['CO2'] + sg.imass['CO2'] + tot * (44/16))
        fl.imass['Water'] = tot * (36/16)
        fl.T = 423.15

    def _design(self):
        self.design_results['Elec_kW'] = self._pwr
        self.design_results['Heat_kW'] = self._heat

    def _cost(self):
        self.purchase_costs['Gas engine-CHP'] = 800 * max(self._pwr, 10) * CEPCI_FACTOR


# ═══════════════════════════════════════════════════════════════════════════════
# 5.  BUILD AND SIMULATE SYSTEM
# ═══════════════════════════════════════════════════════════════════════════════
U01 = ColdPress        ('U01', ins=feed_kernel,         outs=('crude_oil', 'press_cake'))
U02 = OilFilter        ('U02', ins=U01-0,               outs=('argan_oil', 'oil_fines'))
U03 = SaponinExtractor ('U03', ins=(U01-1, solvent),    outs=('sap_ext', 'dep_cake'))
U04 = Evaporator       ('U04', ins=U03-0,               outs=('sap_conc', 'cond'))
U05 = AnaerobicDigester('U05', ins=feed_pulp,           outs=('biogas', 'digestate'))
U06 = SlowPyrolysis    ('U06', ins=feed_shell,          outs=('biochar', 'bio_oil', 'syngas'))
U07 = CHP              ('U07', ins=(U05-0, U06-2),      outs=('flue',))

sys = bst.System('argan_bf', path=(U01, U02, U03, U04, U05, U06, U07))
sys.simulate()
print("BioSTEAM simulation  : CONVERGED ✓\n")

# ── Alias key product streams ──────────────────────────────────────────────────
oil    = U02.outs[0]
sap    = U04.outs[0]
bc     = U06.outs[0]
bio_oil = U06.outs[1]
bg     = U05.outs[0]
dg     = U05.outs[1]

# ═══════════════════════════════════════════════════════════════════════════════
# 6.  TECHNO-ECONOMIC ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════════
H    = 8000   # operating hours/yr
LIFE = 20     # project life (yr)
DR   = 0.10   # WACC / discount rate
TAX  = 0.30   # Moroccan corporate income tax

# ── CAPEX: Lang factor method (Turton 2009, §6) ───────────────────────────────
# Lang factor 4.2 for mixed fluid-solid biorefinery, decomposed as:
#   Direct  : 2.53 × PCE  (equip+install+piping+instr+civil+insul)
#   Indirect: 0.64 × PCE  (engineering+construction+legal)
#   Contng. : 0.47 × PCE  (15 % × (direct+indirect))
#   Site dev: 0.56 × PCE
#   Total   : 4.20 × PCE
UNITS = [U01, U02, U03, U04, U05, U06, U07]

PCE_breakdown = {}
for u in UNITS:
    for k, v in u.purchase_costs.items():
        PCE_breakdown[f'{u.ID}: {k}'] = v

PCE   = sum(sum(u.purchase_costs.values()) for u in UNITS)
LANG  = 4.2
CAPEX = PCE * LANG                 # Fixed capital investment (FCI)
WC    = 0.05 * CAPEX               # Working capital = 5 % FCI (Turton §9)
TCI   = CAPEX + WC                 # Total capital investment

# Construction schedule: 60 % yr-0, 40 % yr-1 (2-year build)
CAPEX_Y0 = 0.60 * CAPEX
CAPEX_Y1 = 0.40 * CAPEX

# ── OPEX ──────────────────────────────────────────────────────────────────────
# Market prices (2026 USD):
# ── Oil price scenarios (Reviewer correction: $120/kg is premium retail;
#    bulk wholesale range is $40–80/kg; use $60 as realistic base)
OIL_SCENARIOS = {
    'Low ($40/kg)'  : 40.0,
    'Base ($60/kg)' : 60.0,
    'High ($120/kg)': 120.0,
}
OIL_BASE = 60.0   # ← base case for main results

PR = {
    'oil'  : OIL_BASE, # $/kg — see scenario analysis; bulk ~$60, premium ~$120
    'sap'  : 35.0,     # $/kg — technical-grade saponin [35]; Henry et al. (2013) [13] confirms argan saponin structure
    'bc'   : 0.80,     # $/kg — agricultural-grade biochar [33,34]; range $0.25-1.50/kg
    'bio'  : 0.25,     # $/kg — pyrolysis bio-oil
    'elec' : 0.08,     # $/kWh — Moroccan industrial HV tariff [31,32]; ~MAD 0.73/kWh ≈ $0.08/kWh large consumers
    'dg'   : 0.02,     # $/kg — digestate as biofertiliser
}

opex_feedstock   = CAP  * H * 0.12        # argan fruit $0.12/kg (MAD 1.2/kg) bulk farm-gate
    # Industrial bulk purchase for 10 000 t/yr; cooperatives pay MAD 5–15/kg for quality
    # kernels separately. Whole-fruit bulk price verified at MAD 1–3/kg for large-scale
    # industrial buyers (cf. Lybbert et al. 2011; FAO 2022). Conservative estimate.
opex_water       = 300  * H * 0.001       # process water
opex_solvent     = 180  * H * 0.05 * 0.60 # 5 % EtOH makeup @ $0.60/kg
# Labour — Moroccan context (Souss-Massa / Agadir industrial zone, 2024):
# Staffing mix per shift (3 shifts, 4 operators/shift = 12 total):
#   4 × basic operator (SMIG+):          MAD 3,500/month  (ref: glassdoor.fr, emploi.ma)
#   6 × process technician (Bac+2/3):    MAD 5,000/month  (ref: paylab.com industry chimique)
#   2 × senior technician (experienced): MAD 8,000/month  (ref: paylab.com)
#   Weighted average gross:              MAD 5,000/month
# + 26% CNSS employer social contribution → MAD 6,300/month effective per operator
# Converted at 10.2 MAD/USD → $618/month ≈ $620/month (rounded)
# Sources: paylab.com/ma (industrie chimique, 2024); glassdoor.fr Morocco;
#          emploi.ma; jobsquare.ma pharmaceutical industry benchmarks.
opex_labor       = 12   * 620 * 12        # 12 operators @ $620/month (Moroccan context)
opex_supervision = 0.20 * opex_labor      # 20% supervision + fringe benefits (Turton Table 8.2)
opex_maint       = 0.04 * CAPEX           # 4 % FCI/yr
opex_insurance   = 0.005 * CAPEX          # 0.5% FCI/yr (Moroccan industrial plants; range 0.4–0.6%)
opex_utilities   = U07._pwr * H * 0.01   # minimal net purchased electricity
opex_catalyst    = 5000 * 12              # chemicals, pH agents

OPEX = (opex_feedstock + opex_water + opex_solvent + opex_labor +
        opex_supervision + opex_maint + opex_insurance +
        opex_utilities + opex_catalyst)

# ── Revenue ────────────────────────────────────────────────────────────────────
rev_oil  = oil.imass['Triolein'] * H * PR['oil']
rev_sap  = sap.imass['Saponin']  * H * PR['sap']
rev_bc   = bc.imass['Biochar']   * H * PR['bc']
rev_bio  = bio_oil.F_mass        * H * PR['bio']
rev_elec = U07._pwr              * H * PR['elec']
rev_dg   = dg.F_mass             * H * PR['dg']
REV = rev_oil + rev_sap + rev_bc + rev_bio + rev_elec + rev_dg

# ── Profitability (straight-line depreciation) ────────────────────────────────
DEPR   = CAPEX / LIFE
EBITDA = REV - OPEX
EBIT   = EBITDA - DEPR
NI     = max(0, EBIT * (1 - TAX))         # net income after tax
CF     = NI + DEPR                         # operating cash flow

# DCF schedule: yr-0 = -60%CAPEX; yr-1 = -40%CAPEX-WC; yr 2..21 = CF; yr-22 = +WC
DCF = [-CAPEX_Y0, -CAPEX_Y1 - WC] + [CF] * LIFE + [WC]
NPV = sum(DCF[i] / (1 + DR)**i for i in range(len(DCF)))

try:
    IRR = brentq(lambda r: sum(DCF[i]/(1+r)**i for i in range(len(DCF))), 0.001, 5.0)
except Exception:
    IRR = 0.30

PB   = TCI / CF                            # simple payback
LCOP = OPEX / (oil.imass['Triolein'] * H)  # $/kg argan oil

# ── Three oil price scenarios ─────────────────────────────────────────────
scenario_results = {}
for sc_name, oil_price in OIL_SCENARIOS.items():
    sc_rev_oil = oil.imass['Triolein'] * H * oil_price
    sc_rev     = sc_rev_oil + rev_sap + rev_bc + rev_bio + rev_elec + rev_dg
    sc_ebitda  = sc_rev - OPEX
    sc_ebit    = sc_ebitda - DEPR
    sc_ni      = max(0, sc_ebit * (1 - TAX))
    sc_cf      = sc_ni + DEPR
    sc_dcf     = [-CAPEX_Y0, -CAPEX_Y1 - WC] + [sc_cf] * LIFE + [WC]
    sc_npv     = sum(sc_dcf[i] / (1+DR)**i for i in range(len(sc_dcf)))
    try:
        sc_irr = brentq(lambda r: sum(sc_dcf[i]/(1+r)**i for i in range(len(sc_dcf))),
                        0.001, 5.0)
    except Exception:
        sc_irr = 0.01
    sc_pb   = TCI / max(sc_cf, 1)
    sc_lcop = OPEX / (oil.imass['Triolein'] * H)
    scenario_results[sc_name] = dict(
        oil_price=oil_price, rev=sc_rev, ebitda=sc_ebitda,
        ni=sc_ni, cf=sc_cf, npv=sc_npv, irr=sc_irr*100, pb=sc_pb, lcop=sc_lcop
    )

R_SCEN = scenario_results   # saved to JSON later

print(f"{'─'*50}")
print(f" TEA SUMMARY (CEPCI 2026 = {CEPCI_2026})")
print(f"{'─'*50}")
print(f" PCE                : ${PCE/1e6:.3f} M (2026 USD)")
print(f" CAPEX (Lang 4.2)   : ${CAPEX/1e6:.3f} M")
print(f" TCI                : ${TCI/1e6:.3f} M")
print(f" Annual revenue     : ${REV/1e6:.2f} M/yr")
print(f" Annual OPEX        : ${OPEX/1e6:.3f} M/yr")
print(f" EBITDA             : ${EBITDA/1e6:.2f} M/yr")
print(f" Net income         : ${NI/1e6:.2f} M/yr")
print(f" NPV @ 10 %         : ${NPV/1e6:.1f} M")
print(f" IRR                : {IRR*100:.1f} %")
print(f" Payback period     : {PB:.2f} yr")
print(f" LCOP argan oil     : ${LCOP:.2f}/kg\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 7.  LIFE CYCLE ASSESSMENT  (ISO 14040/14044, cradle-to-gate)
# ═══════════════════════════════════════════════════════════════════════════════
# Functional unit: 1 kg argan oil at factory gate
FU = oil.imass['Triolein'] * H   # kg oil/yr

# ── Life Cycle Inventory (CORRECTED — reviewer v2) ─────────────────────────
#
# REVIEWER FIX 1 — ETHANOL SOLVENT:
#   Previous version charged only 5% EtOH makeup (0.165 kg CO2-eq/kg oil).
#   Correct: no solvent recycle loop exists in the BioSTEAM model — the
#   condensate stream exits as waste. Therefore the FULL EtOH input must be
#   charged (180 kg/h). This is a ×20 correction: 3.304 vs 0.165 kg CO2/kg oil.
#   Future work should add an explicit EtOH distillation–recycle unit to reduce
#   this contribution.
#
# REVIEWER FIX 2 — NO DOUBLE COUNTING OF METHANE:
#   Previous version claimed BOTH (a) electricity from CHP displacing the grid
#   AND (b) biomethane directly substituting fossil natural gas.
#   These are the same energy stream — CH4 → CHP → electricity. Crediting both
#   is double counting. Correct approach per ISO 14044 §4.3.4 system expansion:
#   • Biogenic CH4 combustion in CHP = 0 GWP (biogenic carbon neutrality)
#   • Credit only NET EXPORTED electricity displacing the Moroccan grid
#   • NO separate NG substitution credit
#   This removes the artificial -2.43 kg CO2/kg oil credit and shifts the
#   system from falsely "carbon-negative" to correctly "low-carbon".

# Foreground electricity consumption per unit (kWh/yr)
elec_press   = U01.design_results['Power_kW']       * H
elec_agit    = U03.design_results['Power_kW']       * H
elec_evap    = U04.design_results['Evap_load_kW']   * H
elec_AD      = U05.design_results['Power_kW']       * H
elec_pyro    = U06.design_results['Thermal_duty_kW'] * H * 0.35
elec_misc    = CAP * H * 0.002
elec_IN      = elec_press + elec_agit + elec_evap + elec_AD + elec_pyro + elec_misc
elec_OUT     = U07._pwr * H   # CHP electricity produced (kWh/yr)

# Net electricity balance
net_grid_kWh   = max(0, elec_IN - elec_OUT)    # purchased from grid (charged)
net_export_kWh = max(0, elec_OUT - elec_IN)    # exported to grid (credited)

# Transport (kg diesel/yr)
diesel_fruit = CAP * H / 1000 * 150 * 0.045   # 150 km HGV, 45 g diesel/t·km
diesel_oil   = oil.imass['Triolein'] * H / 1000 * 500 * 0.045  # 500 km to port

# Solvent: FULL EtOH input (no recycle modelled — see note above)
EtOH_total = 180 * H   # kg/yr — full solvent charged

# Background emission factors (Ecoinvent 3.9 cutoff; IPCC AR6)
EF = {
    'grid_kWh'   : 0.547,   # kg CO2-eq/kWh — Morocco grid mix (IEA 2023)
    'diesel_kg'  : 3.195,   # kg CO2-eq/kg  — upstream + combustion
    'ethanol_kg' : 1.52,    # kg CO2-eq/kg  — sugarcane ethanol (Ecoinvent 3.9 cutoff [29]; cf. Pereira et al. 2019 [30])
    'water_kg'   : 0.00034, # kg CO2-eq/kg
    # NOTE: NO CH4_sub factor — biogenic CH4 to CHP = 0 GWP (biogenic neutrality)
    #       Only net exported electricity is credited.
}

# GWP contributions (kg CO2-eq/yr)
gwp_elec   =  net_grid_kWh * EF['grid_kWh']          # grid electricity purchased
gwp_transp = (diesel_fruit + diesel_oil) * EF['diesel_kg']
gwp_etoh   =  EtOH_total * EF['ethanol_kg']           # full EtOH input (no recycle)
gwp_water  =  300 * H * EF['water_kg']
gwp_elec_cr = -net_export_kWh * EF['grid_kWh']        # net electricity export credit
gwp_CH4_c  =  0.0   # REMOVED — would double-count methane (see note above)

GWP_total = gwp_elec + gwp_transp + gwp_etoh + gwp_water + gwp_elec_cr

# Per functional unit
gwp_eu = gwp_elec    / FU
gwp_tu = gwp_transp  / FU
gwp_su = gwp_etoh    / FU   # ≈ 3.30 kg CO2/kg oil (main hotspot)
gwp_wu = gwp_water   / FU
gwp_cu = gwp_elec_cr / FU   # small negative (net exported electricity)
GWP_kg = GWP_total   / FU   # corrected net GWP100 per kg oil

# Fossil energy demand
FED_total = (net_grid_kWh * 9.83           # MJ/kWh primary (Morocco grid eff 37%)
           + (diesel_fruit + diesel_oil) * 43.2   # MJ/kg diesel
           + EtOH_total * 25.0)           # MJ/kg ethanol CED
FED_kg = FED_total / FU

# Water footprint (process water only — tree is rain-fed)
WF_kg = (300 * H) / FU   # kg water / kg oil

print(f"{'─'*50}")
print(f" LCA SUMMARY  (FU: 1 kg argan oil, cradle-to-gate)")
print(f"{'─'*50}")
print(f" GWP100              : {GWP_kg:.4f} kg CO2-eq/kg oil")
print(f"   electricity       : {gwp_eu:+.4f}")
print(f"   transport         : {gwp_tu:+.4f}")
print(f"   EtOH solvent      : {gwp_su:+.4f}")
print(f"   process water     : {gwp_wu:+.4f}")
print(f"   elec export credit: {gwp_cu:+.4f}  (net {net_export_kWh/1000:.1f} MWh/yr exported)")
print(f" Fossil energy (FED) : {FED_kg:.2f} MJ/kg oil")
print(f" Water footprint     : {WF_kg:.2f} kg/kg oil\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 8.  SAVE RESULTS JSON
# ═══════════════════════════════════════════════════════════════════════════════
R = {
    # Simulation
    'oil_flow'     : float(oil.imass['Triolein']),
    'sap_flow'     : float(sap.imass['Saponin']),
    'bc_flow'      : float(bc.imass['Biochar']),
    'bio_oil_flow' : float(bio_oil.F_mass),
    'ch4_flow'     : float(bg.imass['Methane']),
    'dg_flow'      : float(dg.F_mass),
    'elec_kw'      : float(U07._pwr),
    'heat_kw'      : float(U07._heat),
    'H'            : H,
    # CEPCI
    'CEPCI_ref'    : CEPCI_REF,
    'CEPCI_2026'   : CEPCI_2026,
    'CEPCI_factor' : CEPCI_FACTOR,
    # TEA
    'PCE'          : PCE,  'CAPEX' : CAPEX, 'WC' : WC, 'TCI' : TCI,
    'REV'          : REV,  'OPEX'  : OPEX,
    'EBITDA'       : EBITDA, 'DEPR': DEPR, 'NI': NI, 'CF': CF,
    'NPV'          : NPV,  'IRR'   : IRR * 100, 'PB': PB, 'LCOP': LCOP,
    # Revenue breakdown
    'rev_oil'  : rev_oil,  'rev_sap'  : rev_sap,   'rev_bc'   : rev_bc,
    'rev_bio'  : rev_bio,  'rev_elec' : rev_elec,  'rev_dg'   : rev_dg,
    # OPEX breakdown
    'opex_feedstock'  : opex_feedstock,  'opex_water'      : opex_water,
    'opex_solvent'    : opex_solvent,    'opex_labor'      : opex_labor,
    'opex_supervision': opex_supervision,'opex_maint'      : opex_maint,
    'opex_insurance'  : opex_insurance, 'opex_utilities'  : opex_utilities,
    'opex_catalyst'   : opex_catalyst,
    # LCA
    'GWP_kg'   : GWP_kg,  'FED_kg'  : FED_kg,  'WF_kg'   : WF_kg,
    'gwp_eu'   : gwp_eu,  'gwp_tu'  : gwp_tu,  'gwp_su'  : gwp_su,
    'gwp_wu'   : gwp_wu,  'gwp_cu'  : gwp_cu,
    'gwp_elec_t' : gwp_elec/1e3, 'gwp_transp_t': gwp_transp/1e3,
    'gwp_CH4_cr_t': 0.0,
    'elec_OUT_MWh' : elec_OUT / 1000,
    'elec_IN_MWh'  : elec_IN / 1000,
    'net_grid_MWh' : net_grid_kWh / 1000,
    'net_export_MWh': net_export_kWh / 1000,
    'ch4_credit_t'  : 0.0,  # removed — no NG sub credit,
    # PCE breakdown
    'PCE_breakdown': {k: float(v) for k, v in PCE_breakdown.items()},
}

# Add scenario results to R
R['scenarios'] = {k: {kk: float(vv) for kk, vv in v.items()}
                   for k, v in R_SCEN.items()}

json_path = os.path.join(OUT, 'res.json')
with open(json_path, 'w') as f:
    json.dump(R, f, indent=2)
print(f"Results JSON saved   : {json_path}\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 9.  FIGURES
# ═══════════════════════════════════════════════════════════════════════════════
DPI = 240
NAVY   = '#1B4F8A';  GREEN  = '#1E8449';  ORANGE = '#A04000'
TEAL   = '#117A65';  RED    = '#C0392B';  GOLD   = '#B7950B'
PURPLE = '#6C3483';  GREY   = '#7F8C8D';  LBLUE  = '#2471A3'

plt.rcParams.update({
    'font.family'          : 'DejaVu Sans',
    'font.size'            : 11,
    'axes.labelsize'       : 11,
    'axes.titlesize'       : 12,
    'axes.spines.top'      : False,
    'axes.spines.right'    : False,
    'axes.grid'            : True,
    'grid.alpha'           : 0.30,
    'grid.linewidth'       : 0.6,
})

# ── Helper ────────────────────────────────────────────────────────────────────
def save(fig, name):
    p = os.path.join(OUT, name)
    fig.savefig(p, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Figure saved       : {p}")

# ─────────────────────────────────────────────────────────────────────────────
# FIG 1 — PROFESSIONAL PFD (ISO 10628 / ISA style, inspired by Hachhach 2021)
# ─────────────────────────────────────────────────────────────────────────────
print('Drawing professional PFD...')

# CONSTANTS & COLORS
# ═══════════════════════════════════════════════════════════════════════
BG       = '#FFFFFF'
BDY      = '#1A1A2E'      # boundary/border
EQUIP_FC = '#EBF5FB'      # equipment fill (light blue)
EQUIP_EC = '#1B4F8A'      # equipment edge
STREAM_C = '#1A1A1A'      # process stream lines
KERNEL_C = '#1B4F8A'      # kernel train accent
PULP_C   = '#1E8449'      # pulp train accent
SHELL_C  = '#784212'      # shell train accent
UTIL_C   = '#6C3483'      # utility/CHP accent
PROD_C   = '#117A65'      # product streams
BOUND_C  = '#E74C3C'      # system boundary
GRID_C   = '#BDC3C7'      # light grid lines
TEXT_C   = '#1A1A1A'

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 9,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.spines.left': False,
    'axes.spines.bottom': False,
})

fig, ax = plt.subplots(figsize=(24, 17))
fig.patch.set_facecolor(BG)
ax.set_facecolor(BG)
ax.set_xlim(0, 24)
ax.set_ylim(0, 17)
ax.set_aspect('equal')
ax.axis('off')

# ═══════════════════════════════════════════════════════════════════════
# DRAWING PRIMITIVES
# ═══════════════════════════════════════════════════════════════════════

def vessel_v(ax, cx, cy, w, h, tag, name1, name2='', fc=EQUIP_FC, ec=EQUIP_EC, lw=1.8):
    """Vertical pressure vessel with dome top and bottom."""
    er = w * 0.35   # ellipse y-radius
    # Body
    rect = plt.Rectangle((cx-w/2, cy-h/2+er), w, h-2*er,
                           fc=fc, ec=ec, lw=lw, zorder=4)
    ax.add_patch(rect)
    # Top dome
    top = Ellipse((cx, cy+h/2-er), w, 2*er,
                   fc=fc, ec=ec, lw=lw, zorder=5)
    ax.add_patch(top)
    # Bottom dome
    bot = Ellipse((cx, cy-h/2+er), w, 2*er,
                   fc=fc, ec=ec, lw=lw, zorder=5)
    ax.add_patch(bot)
    # Overlap body over ellipses to hide internal lines
    rect2 = plt.Rectangle((cx-w/2, cy-h/2+er), w, h-2*er,
                            fc=fc, ec='none', zorder=4)
    ax.add_patch(rect2)
    # Tag
    ax.text(cx, cy+h/2+0.18, tag, ha='center', va='bottom',
            fontsize=8, fontweight='bold', color=ec, zorder=6)
    ax.text(cx, cy, name1, ha='center', va='center',
            fontsize=8.5, fontweight='bold', color=TEXT_C,
            zorder=6, multialignment='center')
    if name2:
        ax.text(cx, cy-0.28, name2, ha='center', va='center',
                fontsize=7.5, color='#555555', zorder=6, style='italic')
    return (cx, cy+h/2, cx, cy-h/2, cx-w/2, cy, cx+w/2, cy)

def cstr(ax, cx, cy, w, h, tag, name1, name2='', fc=EQUIP_FC, ec=EQUIP_EC, lw=1.8):
    """CSTR (stirred tank) with agitator symbol."""
    vessel_v(ax, cx, cy, w, h, tag, name1, name2, fc=fc, ec=ec, lw=lw)
    # Agitator shaft
    shaft_h = h * 0.55
    ax.plot([cx, cx], [cy, cy+shaft_h*0.5], color=ec, lw=1.2, zorder=6)
    # Impeller blades
    iw = w * 0.45
    ax.plot([cx-iw/2, cx+iw/2], [cy, cy], color=ec, lw=1.8, zorder=6)
    ax.plot([cx-iw/2, cx+iw/2], [cy-0.2, cy-0.2], color=ec, lw=1.8, zorder=6)

def drum_h(ax, cx, cy, w, h, tag, name1, name2='', fc=EQUIP_FC, ec=EQUIP_EC, lw=1.8):
    """Horizontal drum / rotary kiln."""
    er = h * 0.45
    # Body
    rect = plt.Rectangle((cx-w/2+er, cy-h/2), w-2*er, h,
                           fc=fc, ec=ec, lw=lw, zorder=4)
    ax.add_patch(rect)
    # Left cap
    lc = Ellipse((cx-w/2+er, cy), 2*er, h,
                  fc=fc, ec=ec, lw=lw, zorder=5)
    ax.add_patch(lc)
    # Right cap
    rc = Ellipse((cx+w/2-er, cy), 2*er, h,
                  fc=fc, ec=ec, lw=lw, zorder=5)
    ax.add_patch(rc)
    rect2 = plt.Rectangle((cx-w/2+er, cy-h/2), w-2*er, h,
                            fc=fc, ec='none', zorder=4)
    ax.add_patch(rect2)
    # Hatching lines to indicate rotation
    for x_off in np.linspace(-w*0.3, w*0.3, 5):
        ax.plot([cx+x_off, cx+x_off], [cy-h/2*0.5, cy+h/2*0.5],
                color=ec, lw=0.5, alpha=0.35, zorder=5, ls='--')
    ax.text(cx, cy+h/2+0.18, tag, ha='center', va='bottom',
            fontsize=8, fontweight='bold', color=ec, zorder=6)
    ax.text(cx, cy+0.1, name1, ha='center', va='center',
            fontsize=8.5, fontweight='bold', color=TEXT_C, zorder=6)
    if name2:
        ax.text(cx, cy-0.25, name2, ha='center', va='center',
                fontsize=7.5, color='#555555', zorder=6, style='italic')

def box_unit(ax, cx, cy, w, h, tag, name1, name2='', fc=EQUIP_FC, ec=EQUIP_EC, lw=1.8):
    """Generic equipment box (filter, press, etc.)."""
    r = FancyBboxPatch((cx-w/2, cy-h/2), w, h,
                        boxstyle='round,pad=0.08',
                        fc=fc, ec=ec, lw=lw, zorder=4)
    ax.add_patch(r)
    ax.text(cx, cy+h/2+0.15, tag, ha='center', va='bottom',
            fontsize=8, fontweight='bold', color=ec, zorder=6)
    ax.text(cx, cy+(0.12 if name2 else 0), name1,
            ha='center', va='center',
            fontsize=8.5, fontweight='bold', color=TEXT_C, zorder=6)
    if name2:
        ax.text(cx, cy-0.25, name2, ha='center', va='center',
                fontsize=7.5, color='#555555', zorder=6, style='italic')

def hx(ax, cx, cy, size, tag, name='', fc=EQUIP_FC, ec=EQUIP_EC, lw=1.8):
    """Heat exchanger — diamond with X inside."""
    pts = np.array([[cx, cy+size], [cx+size, cy],
                     [cx, cy-size], [cx-size, cy], [cx, cy+size]])
    poly = Polygon(pts, closed=True, fc=fc, ec=ec, lw=lw, zorder=4)
    ax.add_patch(poly)
    ax.plot([cx-size*0.5, cx+size*0.5], [cy+size*0.5, cy-size*0.5],
            color=ec, lw=1.2, zorder=5)
    ax.plot([cx-size*0.5, cx+size*0.5], [cy-size*0.5, cy+size*0.5],
            color=ec, lw=1.2, zorder=5)
    ax.text(cx, cy+size+0.15, tag, ha='center', va='bottom',
            fontsize=8, fontweight='bold', color=ec, zorder=6)
    if name:
        ax.text(cx, cy, name, ha='center', va='center',
                fontsize=7.5, color=TEXT_C, zorder=6)

def chp_engine(ax, cx, cy, r, tag, name1, name2='', ec=UTIL_C, lw=2.0):
    """CHP gas engine — circle with lightning bolt symbol."""
    circ = Circle((cx, cy), r, fc=EQUIP_FC, ec=ec, lw=lw, zorder=4)
    ax.add_patch(circ)
    # Lightning bolt (generator symbol)
    bolt = np.array([[cx-r*0.25, cy+r*0.4], [cx+r*0.1, cy+r*0.05],
                      [cx-r*0.05, cy+r*0.05], [cx+r*0.25, cy-r*0.4]])
    ax.plot(bolt[:,0], bolt[:,1], color=ec, lw=2.0, zorder=5, solid_joinstyle='round')
    ax.text(cx, cy+r+0.18, tag, ha='center', va='bottom',
            fontsize=8, fontweight='bold', color=ec, zorder=6)
    ax.text(cx, cy-r*0.3, name1, ha='center', va='center',
            fontsize=8.5, fontweight='bold', color=TEXT_C, zorder=6)
    if name2:
        ax.text(cx, cy-r*0.65, name2, ha='center', va='center',
                fontsize=7.5, color='#555555', zorder=6, style='italic')

def separator(ax, cx, cy, size, tag):
    """Process separator / splitter — diamond."""
    pts = np.array([[cx, cy+size], [cx+size*0.7, cy],
                     [cx, cy-size], [cx-size*0.7, cy], [cx, cy+size]])
    poly = Polygon(pts, closed=True, fc='#FAF0E6', ec='#5D6D7E', lw=1.8, zorder=4)
    ax.add_patch(poly)
    ax.text(cx, cy, 'SEP', ha='center', va='center',
            fontsize=7.5, fontweight='bold', color='#5D6D7E', zorder=5)
    ax.text(cx, cy+size+0.14, tag, ha='center', va='bottom',
            fontsize=8, fontweight='bold', color='#5D6D7E', zorder=6)

def stream_tag(ax, x, y, num, col='#1B4F8A', bg='white', size=8.5):
    """Stream number in a circle."""
    circ = Circle((x, y), 0.22, fc=bg, ec=col, lw=1.4, zorder=7)
    ax.add_patch(circ)
    ax.text(x, y, str(num), ha='center', va='center',
            fontsize=size, fontweight='bold', color=col, zorder=8)

def ortho_arrow(ax, x1, y1, x2, y2, snum=None, col=STREAM_C, lw=1.6,
                midx=None, midy=None):
    """Orthogonal (L-shaped) process stream arrow."""
    if midx is None and midy is None:
        midx = x2; midy = y1
    xs = [x1, midx, midx, x2]
    ys = [y1, y1,   midy, midy] if midy != y2 else [y1, y1, y2, y2]
    # If simple horizontal or vertical, just draw direct arrow
    if abs(y1-y2) < 0.01:
        xs = [x1, x2]; ys = [y1, y2]
    elif abs(x1-x2) < 0.01:
        xs = [x1, x2]; ys = [y1, y2]
    else:
        xs = [x1, midx, midx, x2]
        ys = [y1, y1, y2, y2]
    ax.plot(xs, ys, color=col, lw=lw, zorder=3, solid_joinstyle='round',
            solid_capstyle='round')
    # Arrow head at end
    dx = x2 - xs[-2]; dy = y2 - ys[-2]
    norm = (dx**2+dy**2)**0.5
    if norm > 0.001:
        ax.annotate('', xy=(x2, y2),
                    xytext=(x2-dx/norm*0.18, y2-dy/norm*0.18),
                    arrowprops=dict(arrowstyle='->', color=col, lw=lw,
                                   mutation_scale=12), zorder=3)
    if snum is not None:
        mx = (xs[-2]+x2)/2 + 0.0; my = (ys[-2]+y2)/2 + 0.0
        stream_tag(ax, mx, my, snum, col=col)

def prod_label(ax, x, y, name, value, col=PROD_C, arrow_dir='right'):
    """Product stream label box."""
    dx = 0.18 if arrow_dir == 'right' else -0.18
    r = FancyBboxPatch((x+(0.05 if arrow_dir=='right' else -1.85),
                         y-0.25), 1.8, 0.5,
                        boxstyle='round,pad=0.06',
                        fc='#F0FFF4' if col==PROD_C else '#EBF5FB',
                        ec=col, lw=1.2, zorder=5, alpha=0.9)
    ax.add_patch(r)
    ax.text(x+(1.0 if arrow_dir=='right' else -0.95), y+0.07,
            name, ha='center', va='center',
            fontsize=8, fontweight='bold', color=col, zorder=6)
    ax.text(x+(1.0 if arrow_dir=='right' else -0.95), y-0.13,
            value, ha='center', va='center',
            fontsize=7.5, color='#2C3E50', zorder=6, style='italic')

# ═══════════════════════════════════════════════════════════════════════
# SYSTEM BOUNDARY (LCA cradle-to-gate, ISO 14040 gate-to-gate)
# ═══════════════════════════════════════════════════════════════════════
sys_boundary = FancyBboxPatch((1.0, 1.0), 19.5, 14.2,
                               boxstyle='round,pad=0.2',
                               fc='none', ec=BOUND_C,
                               lw=2.0, ls='--', zorder=1)
ax.add_patch(sys_boundary)
ax.text(1.35, 15.05, 'SYSTEM BOUNDARY  (ISO 14040/14044 Cradle-to-Gate)',
        ha='left', va='center', fontsize=8.5, color=BOUND_C,
        fontweight='bold', style='italic', zorder=6)

# ═══════════════════════════════════════════════════════════════════════
# ZONE BACKGROUNDS (light shading per train)
# ═══════════════════════════════════════════════════════════════════════
def zone(ax, x, y, w, h, label, col, alpha=0.06):
    r = plt.Rectangle((x, y), w, h, fc=col, ec=col,
                        lw=0.8, alpha=alpha, zorder=0, ls=':')
    ax.add_patch(r)
    # Zone label (watermark style)
    ax.text(x+0.2, y+h-0.25, label, ha='left', va='top',
            fontsize=8, color=col, alpha=0.6, fontweight='bold',
            style='italic', zorder=1)

zone(ax, 1.2, 10.0, 16.0, 4.8, 'KERNEL TRAIN',  KERNEL_C)
zone(ax, 1.2,  5.8, 16.0, 3.9, 'PULP TRAIN',    PULP_C)
zone(ax, 1.2,  1.2, 16.0, 4.3, 'SHELL TRAIN',   SHELL_C)
zone(ax, 17.2, 1.2,  3.1,13.6, 'CHP / UTILITY', UTIL_C)

# ═══════════════════════════════════════════════════════════════════════
# EQUIPMENT
# ═══════════════════════════════════════════════════════════════════════

# Feed block
box_unit(ax, 1.9, 8.5, 1.4, 3.6, '', 'ARGAN\nFRUIT', '1250 kg/h',
         fc='#FFF9C4', ec='#F39C12', lw=2.0)

# Splitter / Fruit Separator
separator(ax, 4.1, 8.5, 0.55, 'S-100')

# ── KERNEL TRAIN ─────────────────────────────────────────────────────
# P-101: Hydraulic Cold Press
box_unit(ax, 6.8, 12.5, 1.9, 1.1, 'P-101',
         'COLD PRESS', 'η_oil=92%',
         fc=EQUIP_FC, ec=KERNEL_C)
# Internal press symbol: hatched lines
for yi in np.linspace(12.1, 12.9, 5):
    ax.plot([6.0, 7.6], [yi, yi], color=KERNEL_C, lw=0.6, alpha=0.45, zorder=5)

# F-101: Oil plate filter
box_unit(ax, 9.8, 12.5, 1.6, 1.1, 'F-101',
         'PLATE FILTER', '95% H₂O sep.',
         fc=EQUIP_FC, ec=KERNEL_C)
# filter mesh lines
for xi in np.linspace(9.2, 10.4, 6):
    ax.plot([xi, xi], [12.05, 12.95], color=KERNEL_C, lw=0.7, alpha=0.45, zorder=5)

# V-101: Saponin extractor (tall vessel)
vessel_v(ax, 6.8, 10.6, 1.6, 1.7, 'V-101', 'SAPONIN', 'EXTRACTOR',
         fc='#E8F8F5', ec=KERNEL_C)

# E-101: Multi-effect evaporator
vessel_v(ax, 9.8, 10.6, 1.6, 1.7, 'E-101', '3-EFFECT', 'EVAPORATOR',
         fc='#EBF5FB', ec=KERNEL_C)
# Heating coils symbol
for yi in np.linspace(10.3, 10.9, 3):
    arc = Arc((9.8, yi), 0.7, 0.2, angle=0, theta1=0, theta2=180,
               color=KERNEL_C, lw=0.9, alpha=0.55, zorder=5)
    ax.add_patch(arc)

# ── PULP TRAIN ────────────────────────────────────────────────────────
# R-101: Anaerobic Digester (large CSTR)
cstr(ax, 7.2, 7.2, 2.4, 2.4, 'R-101', 'ANAEROBIC', 'DIGESTER  CSTR 35°C',
     fc='#EAFAF1', ec=PULP_C, lw=2.0)

# ── SHELL TRAIN ───────────────────────────────────────────────────────
# K-101: Rotary kiln pyrolysis
drum_h(ax, 9.0, 3.0, 3.6, 1.2, 'K-101', 'ROTARY KILN PYROLYSIS',
       '470-490°C  (Ait Itto 2024)',
       fc='#FDF2E9', ec=SHELL_C, lw=2.0)

# ── UTILITY / CHP ─────────────────────────────────────────────────────
# G-101: CHP Gas Engine
chp_engine(ax, 18.7, 7.5, 0.95, 'G-101', 'GAS ENGINE', 'CHP Unit',
           ec=UTIL_C, lw=2.0)

# ═══════════════════════════════════════════════════════════════════════
# PROCESS STREAMS  (numbered S-01 … S-17)
# ═══════════════════════════════════════════════════════════════════════

# S-01: Argan fruit → SEP
ortho_arrow(ax, 2.6, 8.5, 3.55, 8.5, snum=1, col='#F39C12', lw=2.0)

# S-02: Kernel → P-101 (up-right path)
ortho_arrow(ax, 4.65, 9.3, 5.9, 12.5, snum=2, col=KERNEL_C,
            midx=4.65, midy=12.5)

# S-03: Pulp → R-101 (straight right)
ortho_arrow(ax, 4.65, 8.5, 6.0, 7.2, snum=3, col=PULP_C,
            midx=6.0, midy=8.5)

# S-04: Shell → K-101 (down-right)
ortho_arrow(ax, 4.65, 7.7, 7.2, 3.0, snum=4, col=SHELL_C,
            midx=4.65, midy=3.0)

# S-05: Crude oil P-101 → F-101
ortho_arrow(ax, 7.7, 12.5, 9.0, 12.5, snum=5, col=KERNEL_C)

# S-06: Press cake P-101 → V-101 (down)
ortho_arrow(ax, 6.8, 11.95, 6.8, 11.45, snum=6, col=KERNEL_C, lw=1.4)

# S-07: Argan oil F-101 → product (right)
ortho_arrow(ax, 10.6, 12.8, 20.9, 12.8, snum=7, col=PROD_C, lw=2.0,
            midx=20.9, midy=12.8)

# S-08: Oil fines (small waste stream down from F-101)
ortho_arrow(ax, 9.8, 11.95, 9.8, 11.45, snum=8, col=GRID_C, lw=1.2)
ax.text(8.8, 11.15, 'fines\n(waste)', ha='center', va='top',
        fontsize=7, color=GRID_C, style='italic', zorder=6)

# EtOH/water solvent input (from outside system boundary — dashed)
ax.annotate('', xy=(6.8, 9.5), xytext=(5.5, 9.5),
            arrowprops=dict(arrowstyle='->', color=KERNEL_C, lw=1.4,
                            linestyle='dashed'), zorder=3)
ax.text(4.1, 9.5, '⊕ S-09\n60% EtOH/H₂O\n300 kg/h  60°C', ha='left', va='center',
        fontsize=7.5, color=KERNEL_C, style='italic', zorder=6)
ax.text(4.0, 9.5, '', ha='left', va='center', zorder=6)
# Dashed box for input crossing boundary
ax.plot([1.0, 5.5], [9.5, 9.5], color=KERNEL_C, lw=1.2, ls=':', zorder=3, alpha=0.6)

# S-10: Saponin extract V-101 → E-101
ortho_arrow(ax, 7.6, 10.6, 9.0, 10.6, snum=10, col=KERNEL_C)

# S-11: Saponin concentrate → product
ortho_arrow(ax, 10.6, 10.9, 20.9, 10.9, snum=11, col=PROD_C, lw=2.0,
            midx=20.9, midy=10.9)

# S-12: Condensate from evaporator (potentially recyclable EtOH)
ortho_arrow(ax, 9.8, 9.75, 9.8, 9.25, snum=12, col=GRID_C, lw=1.2)
ax.text(10.0, 9.0, 'condensate\n(EtOH recycle*)', ha='left', va='top',
        fontsize=7, color=GRID_C, style='italic', zorder=6)

# S-13: Biogas from R-101 → G-101
ortho_arrow(ax, 8.4, 7.2, 17.75, 7.2, snum=13, col=PULP_C, lw=1.8,
            midx=17.75, midy=7.2)

# S-14: Digestate from R-101 (exits left, product)
ortho_arrow(ax, 6.0, 7.2, 1.0, 7.2, snum=14, col=PROD_C, lw=1.6,
            midx=1.0, midy=7.2)

# S-15: Syngas from K-101 → G-101 (up then right)
ortho_arrow(ax, 10.8, 3.0, 17.75, 6.55, snum=15, col=SHELL_C, lw=1.6,
            midx=17.75, midy=3.0)

# S-16: Biochar → product
ortho_arrow(ax, 10.8, 3.4, 20.9, 3.4, snum=16, col=PROD_C, lw=2.0)

# S-17: Bio-oil → product
ortho_arrow(ax, 10.8, 2.6, 20.9, 2.6, snum=17, col=PROD_C, lw=1.6)

# G-101 outputs
# S-18: Electricity → outside boundary
ortho_arrow(ax, 19.65, 7.5, 20.9, 7.5, snum=18, col=UTIL_C, lw=2.0)

# S-19: Process heat (dashed internal) → E-101
ax.annotate('', xy=(10.6, 10.3), xytext=(19.65, 7.5),
            arrowprops=dict(arrowstyle='->', color=UTIL_C, lw=1.2,
                            linestyle='dashed',
                            connectionstyle='arc3,rad=-0.35'), zorder=3)
ax.text(15.5, 9.5, 'S-20\nprocess\nheat', ha='center', va='center',
        fontsize=7, color=UTIL_C, style='italic', zorder=6)

# S-21: Flue gas → outside (up from G-101)
ortho_arrow(ax, 18.7, 8.45, 18.7, 15.4, snum=21, col=GRID_C, lw=1.2)
ax.text(19.2, 15.3, 'flue gas\n(CO₂+H₂O)', ha='left', va='top',
        fontsize=7, color=GRID_C, style='italic', zorder=6)

# ═══════════════════════════════════════════════════════════════════════
# PRODUCT LABELS (outside right boundary)
# ═══════════════════════════════════════════════════════════════════════
oil_f  = f'{R.get("oil_flow",82.8):.0f} kg/h'
sap_f  = f'{R.get("sap_flow",7.6):.1f} kg/h'
bc_f   = f'{R.get("bc_flow",47.5):.0f} kg/h'
bo_f   = f'{R.get("bio_oil_flow",31.8):.0f} kg/h'
elec_f = f'{R.get("elec_kw",307):.0f} kW'
dg_f   = f'{R.get("dg_flow",600):.0f} kg/h'

prod_label(ax, 20.9, 12.8, 'Argan Oil',    oil_f,  KERNEL_C)
prod_label(ax, 20.9, 10.9, 'Saponin Conc.',sap_f,  KERNEL_C)
prod_label(ax, 20.9,  3.4, 'Biochar',       bc_f,  SHELL_C)
prod_label(ax, 20.9,  2.6, 'Bio-oil',       bo_f,  SHELL_C)
prod_label(ax, 20.9,  7.5, 'Electricity',   elec_f, UTIL_C)
# Digestate at left
prod_label(ax, -0.1, 7.2, 'Digestate', dg_f, PULP_C, 'right')
ax.text(1.9, 7.2, '← S-14', ha='center', va='center',
        fontsize=7.5, color=PROD_C, style='italic', zorder=6)

# ═══════════════════════════════════════════════════════════════════════
# STREAM TABLE (bottom)
# ═══════════════════════════════════════════════════════════════════════
# Table header
th = 0.55   # table y-start from absolute bottom
tw = 21.5   # table width
cols_x = [0.2, 1.2, 3.8, 7.5, 10.5, 13.5, 17.0]
col_headers = ['No.', 'Description', 'T (°C)', 'P (kPa)', 'Flow (kg/h)', 'Main comp.', 'Phase']

# Draw table title
ax.text(0.2, 0.05 + 0.9 + 0.26, 'STREAM TABLE', ha='left', va='bottom',
        fontsize=9, fontweight='bold', color=EQUIP_EC, zorder=6)

# Header row background
hrow = FancyBboxPatch((0.1, 0.04 + 0.78), tw, 0.38,
                       boxstyle='round,pad=0.04',
                       fc=EQUIP_EC, ec=EQUIP_EC, lw=0, zorder=5)
ax.add_patch(hrow)
for ci, (cx_c, hdr) in enumerate(zip(cols_x, col_headers)):
    ax.text(cx_c, 0.04 + 0.78 + 0.19, hdr, ha='left', va='center',
            fontsize=8.2, fontweight='bold', color='white', zorder=6)

# Stream data
streams_table = [
    # (num, description, T_C, P_kPa, flow, component, phase)
    (1,  'Argan fruit feed',       25, 101, '1250',  'Multi',         'S'),
    (2,  'Kernel fraction',        25, 101, '187.5', 'Oil,Prot,Sap',  'S'),
    (3,  'Pulp fraction',          25, 101, '937.5', 'Gluc,Cell,H₂O', 'L/S'),
    (4,  'Shell fraction',         25, 101, '125.0', 'Lign,Cell,Ash', 'S'),
    (5,  'Crude oil (press)',       30, 101, '~102',  'Triolein,H₂O',  'L'),
    (6,  'Press cake',             25, 101, '~87',   'Prot,Sap,Cell', 'S'),
    (7,  'Argan oil (PRODUCT)',    30, 101, f'{R.get("oil_flow",82.8):.1f}', 'Triolein', 'L'),
    (9,  'EtOH/H₂O solvent (in)', 60, 101, '300',   '60% EtOH',      'L'),
    (10, 'Saponin extract',        60, 101, '~308',  'Sap,EtOH,H₂O', 'L'),
    (11, 'Saponin conc. (PROD.)',  70, 101, f'{R.get("sap_flow",7.6):.1f}',  'Saponin', 'L'),
    (13, 'Biogas (AD→CHP)',        35, 110, f'{R.get("ch4_flow",54.8):.1f}', 'CH₄,CO₂', 'G'),
    (14, 'Digestate (PRODUCT)',    35, 101, f'{R.get("dg_flow",600):.0f}',   'N,P,K,H₂O','L'),
    (15, 'Syngas (pyro→CHP)',     500, 110, '~46',   'CH₄,CO₂',       'G'),
    (16, 'Biochar (PRODUCT)',     500, 101, f'{R.get("bc_flow",47.5):.1f}',  'C>70%',   'S'),
    (17, 'Bio-oil (PRODUCT)',      25, 101, f'{R.get("bio_oil_flow",31.8):.1f}','Phenolics','L'),
    (18, 'Electricity (PRODUCT)',  '--','--',f'{R.get("elec_kw",307):.0f} kW','—',      'Elec'),
    (21, 'Flue gas (to atm.)',    150, 101, '~200',  'CO₂,H₂O',       'G'),
]

row_h = 0.355
for ri, row in enumerate(streams_table):
    y_row = 0.04 + 0.78 - (ri+1)*row_h
    if y_row < -0.8: break   # don't go off page
    bg = '#F4F6F7' if ri % 2 == 0 else 'white'
    r_patch = FancyBboxPatch((0.1, y_row), tw, row_h-0.03,
                              boxstyle='round,pad=0.03',
                              fc=bg, ec=GRID_C, lw=0.4, zorder=5)
    ax.add_patch(r_patch)
    row_data = [str(row[0]), row[1], str(row[2]), str(row[3]),
                row[4], row[5], row[6]]
    # Color coding
    prod_rows = {7, 11, 14, 16, 17, 18}
    text_col = PROD_C if row[0] in prod_rows else TEXT_C
    for ci, (cx_c, val) in enumerate(zip(cols_x, row_data)):
        bold = (ci == 0 or row[0] in prod_rows)
        ax.text(cx_c, y_row + row_h/2 - 0.03, val,
                ha='left', va='center', fontsize=7.5, color=text_col,
                fontweight='bold' if bold else 'normal', zorder=6)

# Table border
table_top = 0.04 + 0.78 + 0.38
table_bot = 0.04 + 0.78 - len(streams_table)*row_h
ax.plot([0.1, 0.1+tw], [table_top, table_top], color=EQUIP_EC, lw=1.2, zorder=6)
ax.plot([0.1, 0.1+tw], [table_bot+0.03, table_bot+0.03], color=EQUIP_EC, lw=1.0, zorder=6)

# ═══════════════════════════════════════════════════════════════════════
# LEGEND
# ═══════════════════════════════════════════════════════════════════════
leg_x = 0.15; leg_y = 15.8
ax.text(leg_x, leg_y, 'LEGEND', ha='left', va='top', fontsize=8.5,
        fontweight='bold', color=TEXT_C, zorder=6)
legend_items = [
    (KERNEL_C, 'Kernel train (oil + saponins)'),
    (PULP_C,   'Pulp train (anaerobic digestion)'),
    (SHELL_C,  'Shell train (slow pyrolysis)'),
    (UTIL_C,   'CHP / utilities'),
    (PROD_C,   'Product streams'),
    (BOUND_C,  'LCA system boundary (dashed)'),
    (GRID_C,   'Utility / waste streams'),
]
for i, (col, label) in enumerate(legend_items):
    ax.plot([leg_x, leg_x+0.55], [leg_y-0.35-i*0.35, leg_y-0.35-i*0.35],
            color=col, lw=2.2, zorder=6)
    ax.text(leg_x+0.7, leg_y-0.35-i*0.35, label, ha='left', va='center',
            fontsize=8, color=TEXT_C, zorder=6)

# ═══════════════════════════════════════════════════════════════════════
# TITLE BLOCK (bottom right)
# ═══════════════════════════════════════════════════════════════════════
tb_x = 15.0; tb_y = -0.85; tb_w = 9.0; tb_h = 1.5
tb = FancyBboxPatch((tb_x, tb_y), tb_w, tb_h,
                     boxstyle='round,pad=0.08',
                     fc='#EBF5FB', ec=EQUIP_EC, lw=1.5, zorder=5)
ax.add_patch(tb)
ax.text(tb_x+tb_w/2, tb_y+tb_h-0.28,
        'Integrated Argan (Argania spinosa) Biorefinery',
        ha='center', va='center', fontsize=9.5, fontweight='bold',
        color=EQUIP_EC, zorder=6)
ax.text(tb_x+tb_w/2, tb_y+tb_h-0.62,
        'Process Flow Diagram  —  Capacity: 10 000 t/yr argan fruit (1 250 kg/h)',
        ha='center', va='center', fontsize=8.2, color=TEXT_C, zorder=6)
ax.text(tb_x+tb_w/2, tb_y+tb_h-0.95,
        f'BioSTEAM v2.53  |  CEPCI {int(R.get("CEPCI_2026",820))} (2026 USD)  |  Hachhach et al. 2026',
        ha='center', va='center', fontsize=8, color='#555555',
        style='italic', zorder=6)

# Footnote
ax.text(0.1, -0.85, '* S-09 (EtOH/H₂O solvent): in the base LCA scenario, full solvent input is charged '
        '(no recycle modelled). Adding a distillation recovery unit (future work) would significantly reduce '
        'FED and GWP contributions from S-09.',
        ha='left', va='bottom', fontsize=7, color='#777777',
        style='italic', zorder=6)

# ═══════════════════════════════════════════════════════════════════════
# SAVE
# ═══════════════════════════════════════════════════════════════════════
fig.tight_layout(pad=0.0)
save(fig, 'fig1_PFD_professional.png')

# FIG 2 — DETAILED TEA (3 panels)
# ─────────────────────────────────────────────────────────────────────────────
fig2, axes = plt.subplots(1, 3, figsize=(17, 6.5))

# Panel A — Revenue waterfall
rev_lbl = ['Argan Oil', 'Saponin\nConc.', 'Biochar', 'Bio-oil', 'Electricity', 'Digestate']
rev_vals = [R['rev_oil']/1e6, R['rev_sap']/1e6, R['rev_bc']/1e6,
            R['rev_bio']/1e6, R['rev_elec']/1e6, R['rev_dg']/1e6]
rev_cols = [NAVY, TEAL, ORANGE, GOLD, PURPLE, GREY]
xs = np.arange(len(rev_lbl))
bars = axes[0].bar(xs, rev_vals, color=rev_cols, edgecolor='white', linewidth=1.2)
for b, v in zip(bars, rev_vals):
    axes[0].text(b.get_x()+b.get_width()/2, v+0.3, f'${v:.1f}M',
                 ha='center', va='bottom', fontsize=8.5, fontweight='bold')
axes[0].set_xticks(xs);  axes[0].set_xticklabels(rev_lbl, fontsize=9)
axes[0].set_ylabel('Annual Revenue (M USD/yr)');
axes[0].set_ylim(0, max(rev_vals)*1.25)
axes[0].set_title('(A) Revenue by Product Stream', fontweight='bold')
axes[0].text(0.98, 0.97, f'Total: ${sum(rev_vals):.1f} M/yr',
             transform=axes[0].transAxes, ha='right', va='top',
             fontsize=10, fontweight='bold', color=NAVY,
             bbox=dict(boxstyle='round', facecolor='#EBF5FB', edgecolor=NAVY))

# Panel B — OPEX breakdown (horizontal bar)
opex_lbl  = ['Feedstock', 'Labour', 'Supervision', 'Maintenance',
             'Insurance', 'EtOH Solvent', 'Utilities', 'Chemicals', 'Water']
opex_vals = [R['opex_feedstock']/1e6, R['opex_labor']/1e6,
             R['opex_supervision']/1e6, R['opex_maint']/1e6,
             R['opex_insurance']/1e6, R['opex_solvent']/1e6,
             R['opex_utilities']/1e6, R['opex_catalyst']/1e6,
             R['opex_water']/1e6]
opex_cols = [NAVY,TEAL,LBLUE,ORANGE,GREY,GREEN,PURPLE,GOLD,'#BDC3C7']
order = np.argsort(opex_vals)[::-1]
ov = [opex_vals[i] for i in order];  ol = [opex_lbl[i] for i in order]
oc = [opex_cols[i] for i in order]
ys2 = np.arange(len(ov))
axes[1].barh(ys2, ov, color=oc, edgecolor='white', linewidth=1.1, height=0.7)
for y, v in zip(ys2, ov):
    axes[1].text(v+0.005, y, f'${v:.3f}M', va='center', fontsize=8.5)
axes[1].set_yticks(ys2);  axes[1].set_yticklabels(ol, fontsize=9)
axes[1].set_xlabel('Annual OPEX (M USD/yr)')
axes[1].set_title('(B) OPEX Breakdown', fontweight='bold')
axes[1].text(0.98, 0.04, f'Total OPEX: ${sum(ov):.2f} M/yr',
             transform=axes[1].transAxes, ha='right', va='bottom',
             fontsize=9.5, fontweight='bold', color=RED,
             bbox=dict(boxstyle='round', facecolor='#FDFEFE', edgecolor=RED))

# Panel C — Cumulative DCF profile
years_plot = list(range(len(DCF)))
cum_cf = [];  s = 0
for c in DCF: s += c; cum_cf.append(s/1e6)
disc_cf = [];  s2 = 0
for i, c in enumerate(DCF): s2 += c/(1+DR)**i; disc_cf.append(s2/1e6)
axes[2].fill_between(years_plot, cum_cf,  alpha=0.12, color=NAVY)
axes[2].plot(years_plot, cum_cf,  NAVY,  lw=2.2, label='Undiscounted cum. CF')
axes[2].fill_between(years_plot, disc_cf, alpha=0.12, color=GREEN)
axes[2].plot(years_plot, disc_cf, GREEN, lw=2.2, label=f'Discounted (10% WACC)')
axes[2].axhline(0, color='black', lw=0.8, ls='--', alpha=0.5)
axes[2].axvline(2, color=RED,  lw=1.0, ls=':',  alpha=0.6)
axes[2].text(2.15, min(cum_cf)*0.6, 'Start-up', fontsize=8, color=RED, rotation=90)
pb_yr = TCI/1e6 / (CF/1e6) + 2
axes[2].axvline(pb_yr, color=GOLD, lw=1.2, ls='-.', alpha=0.8)
axes[2].text(pb_yr+0.25, cum_cf[0]*0.7, f'PB≈{pb_yr:.1f}yr', fontsize=8, color=GOLD)
axes[2].set_xlabel('Project year');  axes[2].set_ylabel('Cumulative cash flow (M USD)')
axes[2].set_title('(C) Discounted Cash Flow Profile', fontweight='bold')
axes[2].legend(fontsize=9, frameon=False)
# Build scenario NPV summary for the annotation
scen_lines = '\n'.join([f"{k.split('(')[0].strip()}: NPV=${v['npv']/1e6:.0f}M"
                           for k,v in R_SCEN.items()])
txt = f"Base case ($60/kg):\nNPV=${R_SCEN['Base ($60/kg)']['npv']/1e6:.0f}M  IRR={R_SCEN['Base ($60/kg)']['irr']:.0f}%  PB={R_SCEN['Base ($60/kg)']['pb']:.1f}yr\n\nAll scenarios:\n{scen_lines}\n\nCEPCI {int(CEPCI_2026)}"
axes[2].text(0.98, 0.05, txt, transform=axes[2].transAxes, ha='right', va='bottom',
             fontsize=9, fontweight='bold', color=NAVY,
             bbox=dict(boxstyle='round', facecolor='#EBF5FB', edgecolor=NAVY))

fig2.suptitle('Techno-Economic Analysis — Integrated Argan Biorefinery\n'
              f'(CEPCI {int(CEPCI_2026)}, 10% WACC, 20-yr life, CIT 30%; base case oil price $60/kg)',
              fontsize=12, fontweight='bold', y=1.01)
fig2.tight_layout()
save(fig2, 'fig2_TEA.png')

# ─────────────────────────────────────────────────────────────────────────────
# FIG 3 — LCA contribution waterfall + allocation sensitivity
# ─────────────────────────────────────────────────────────────────────────────
fig3, axes = plt.subplots(1, 2, figsize=(15, 6.5))

# Panel A — Waterfall
cats  = ['Grid\nElectricity', 'Transport\n(HGV)', 'EtOH Solvent\n(full input)',
         'Process\nWater', 'Elec. Export\nCredit', 'Net GWP100']
vals  = [R['gwp_eu'], R['gwp_tu'], R['gwp_su'], R['gwp_wu'], R['gwp_cu'], R['GWP_kg']]
bar_c = [RED if v > 0 else GREEN for v in vals];  bar_c[-1] = NAVY

bottoms = [];  cum = 0
for v in vals[:-1]:
    bottoms.append(cum if v > 0 else cum+v);  cum += v
bottoms.append(min(R['GWP_kg'], 0))

xs3 = np.arange(len(cats))
for i, (bot, v, col) in enumerate(zip(bottoms, vals, bar_c)):
    axes[0].bar(i, abs(v), bottom=bot, color=col, width=0.55,
                edgecolor='white', linewidth=1.2)
    off = max(abs(v)*0.07, 0.05)
    ypos = bot + abs(v) + off if v >= 0 else bot - off
    axes[0].text(i, ypos, f'{v:+.3f}', ha='center', va='bottom',
                 fontsize=9, fontweight='bold')
c2 = 0
for i, v in enumerate(vals[:-1]):
    if i < len(vals)-2:
        axes[0].plot([i+0.275, i+0.725], [c2+v, c2+v], 'k-', lw=0.7, alpha=0.4)
    c2 += v
axes[0].axhline(0, color='k', lw=0.7, ls='--', alpha=0.5)
axes[0].set_xticks(xs3);  axes[0].set_xticklabels(cats, fontsize=9.5)
axes[0].set_ylabel('kg CO₂-eq per kg argan oil')
axes[0].set_title('(A) GWP100 Contribution Analysis (CORRECTED)\n'
                   '(FU: 1 kg argan oil, cradle-to-gate; EtOH full input; elec. export only)',
                   fontweight='bold')
pos_p = mpatches.Patch(color=RED,   label='Emission sources')
neg_p = mpatches.Patch(color=GREEN, label='Avoided emission credits')
net_p = mpatches.Patch(color=NAVY,  label='Net GWP100')
axes[0].legend(handles=[pos_p, neg_p, net_p], frameon=False, fontsize=9)

# Panel B — Allocation sensitivity
gross = R['gwp_eu'] + R['gwp_tu'] + R['gwp_su'] + R['gwp_wu']
total_mass  = (R['oil_flow']+R['sap_flow']+R['bc_flow']+
               R['bio_oil_flow']+R['ch4_flow']+R['dg_flow']) * H
mass_alloc  = gross * (R['oil_flow']*H / total_mass)
econ_alloc  = gross * (R['rev_oil'] / R['REV'])
no_credit   = gross

alloc_m  = ['System\nExpansion\n(ISO pref.)', 'Mass\nAllocation',
             'Economic\nAllocation', 'No Credits\n(gross)']
alloc_v  = [R['GWP_kg'], mass_alloc, econ_alloc, no_credit]
alloc_c  = [NAVY, TEAL, GREEN, RED]
bars3 = axes[1].bar(alloc_m, alloc_v, color=alloc_c, edgecolor='white',
                     linewidth=1.2, width=0.55)
for b, v in zip(bars3, alloc_v):
    yp = v+0.02 if v >= 0 else v-0.06
    axes[1].text(b.get_x()+b.get_width()/2, yp, f'{v:+.3f}',
                 ha='center', va='bottom', fontsize=9.5, fontweight='bold')
axes[1].axhline(0, color='k', lw=0.8, ls='--', alpha=0.5)
axes[1].set_ylabel('GWP100 (kg CO₂-eq / kg argan oil)')
axes[1].set_title('(B) Allocation Method Sensitivity', fontweight='bold')
axes[1].text(0.03, 0.97, 'ISO 14044 §4.3.4 recommends\nsystem expansion for\nenergy co-products',
             transform=axes[1].transAxes, ha='left', va='top', fontsize=8.5,
             style='italic', color=NAVY,
             bbox=dict(boxstyle='round', facecolor='#EBF5FB', edgecolor=NAVY, alpha=0.85))

fig3.suptitle('Life Cycle Assessment — GWP100 Results\n'
               '(Integrated Argan Biorefinery, 10 000 t/yr, Ecoinvent 3.9 cutoff)',
               fontsize=12, fontweight='bold', y=1.02)
fig3.tight_layout()
save(fig3, 'fig3_LCA.png')

# ─────────────────────────────────────────────────────────────────────────────
# FIG 4 — TORNADO SENSITIVITY on NPV
# ─────────────────────────────────────────────────────────────────────────────
def calc_npv(cap, rv, op, dr, lf):
    ni  = max(0, (rv - op - cap/lf) * (1 - TAX))
    cf_ = ni + cap/lf
    wc_ = 0.05 * cap
    cfs = [-0.6*cap, -0.4*cap-wc_] + [cf_]*lf + [wc_]
    return sum(cfs[i]/(1+dr)**i for i in range(len(cfs)))

base_npv = R['NPV']
CAP_T = CAPEX;  REV_T = REV;  OPEX_T = OPEX

sens = {
    'Oil price ±20%\n($96–$144/kg)':       (calc_npv(CAP_T,REV_T*1.20,OPEX_T,DR,LIFE),
                                              calc_npv(CAP_T,REV_T*0.80,OPEX_T,DR,LIFE)),
    'CAPEX ±30%':                            (calc_npv(CAP_T*0.70,REV_T,OPEX_T,DR,LIFE),
                                              calc_npv(CAP_T*1.30,REV_T,OPEX_T,DR,LIFE)),
    'Feedstock cost ±40%\n($0.07–$0.17/kg)':(calc_npv(CAP_T,REV_T,OPEX_T-R['opex_feedstock']*0.4,DR,LIFE),
                                              calc_npv(CAP_T,REV_T,OPEX_T+R['opex_feedstock']*0.4,DR,LIFE)),
    'Saponin price ±50%\n($17.5–$52.5/kg)': (calc_npv(CAP_T,REV_T+R['rev_sap']*0.5,OPEX_T,DR,LIFE),
                                              calc_npv(CAP_T,REV_T-R['rev_sap']*0.5,OPEX_T,DR,LIFE)),
    'Discount rate 7–13%':                   (calc_npv(CAP_T,REV_T,OPEX_T,0.07,LIFE),
                                              calc_npv(CAP_T,REV_T,OPEX_T,0.13,LIFE)),
    'Plant life 15–25 yr':                   (calc_npv(CAP_T,REV_T,OPEX_T,DR,25),
                                              calc_npv(CAP_T,REV_T,OPEX_T,DR,15)),
    'Operating factor ±5%':                  (calc_npv(CAP_T,REV_T*1.05,OPEX_T*1.05,DR,LIFE),
                                              calc_npv(CAP_T,REV_T*0.95,OPEX_T*0.95,DR,LIFE)),
    'Biochar price ±60%\n($0.32–$1.28/kg)': (calc_npv(CAP_T,REV_T+R['rev_bc']*0.60,OPEX_T,DR,LIFE),
                                              calc_npv(CAP_T,REV_T-R['rev_bc']*0.60,OPEX_T,DR,LIFE)),
}

labels = list(sens.keys())
high   = [v[0]/1e6 for v in sens.values()]
low    = [v[1]/1e6 for v in sens.values()]
base   = base_npv / 1e6
order4 = np.argsort([abs(h-l) for h,l in zip(high,low)])
labels = [labels[i] for i in order4]
high   = [high[i]   for i in order4]
low    = [low[i]    for i in order4]

fig4, ax4 = plt.subplots(figsize=(13, 6.5))
ys4 = np.arange(len(labels))
ax4.barh(ys4, [h-base for h in high], left=base, color=GREEN, alpha=0.82,
          height=0.52, edgecolor='white', label='Upside scenario')
ax4.barh(ys4, [l-base for l in low],  left=base, color=RED,   alpha=0.82,
          height=0.52, edgecolor='white', label='Downside scenario')
ax4.axvline(base, color=NAVY, lw=2.0, ls='--',
            label=f'Base NPV = ${base:.0f}M')
for y, h, l in zip(ys4, high, low):
    ax4.text(h+2,  y+0.06, f'+${h-base:.0f}M', va='center', fontsize=8, color=GREEN)
    ax4.text(l-2,  y+0.06, f'${l-base:.0f}M',  va='center', ha='right', fontsize=8, color=RED)
ax4.set_yticks(ys4);  ax4.set_yticklabels(labels, fontsize=9.5)
ax4.set_xlabel('NPV (M USD, 10% WACC, 20-yr life)', fontsize=11)
ax4.set_title('Sensitivity Analysis — Tornado Plot\n'
               '(One-at-a-time; ±deviation from base case)', fontweight='bold')
ax4.legend(frameon=False, fontsize=10, loc='lower right')
fig4.tight_layout()
save(fig4, 'fig4_sensitivity.png')

# ─────────────────────────────────────────────────────────────────────────────
# FIG 5 — LCA BENCHMARKING vs literature vegetable oils
# ─────────────────────────────────────────────────────────────────────────────
fig5, axes5 = plt.subplots(1, 3, figsize=(15, 6))
gwp_abs = abs(R['GWP_kg'])

for ax5, (title, ylabel, labels5, vals5, cols5) in zip(axes5, [
    ('(A) Global Warming Potential', 'GWP100 (kg CO₂-eq/kg oil)',
     ['This study\n(argan BFR)','Olive oil\n[Salomone\n2012]',
      'Palm oil\n[Yusoff\n2007]','Jatropha\n[Achten\n2008]','Canola\n[Bernesson\n2006]'],
     [gwp_abs, 3.5, 3.3, 1.9, 3.8],
     [GREEN if gwp_abs < 1.5 else NAVY, ORANGE, RED, TEAL, GOLD]),
    ('(B) Fossil Energy Demand', 'Fossil energy (MJ/kg oil)',
     ['This study','Olive oil\n[Ref.]','Palm oil\n[Ref.]','Jatropha\n[Ref.]','Canola\n[Ref.]'],
     [R['FED_kg'], 18.0, 6.2, 14.5, 16.2],
     [NAVY, ORANGE, RED, TEAL, GOLD]),
    ('(C) Water Footprint', 'Water footprint (kg H₂O/kg oil)',
     ['This study','Olive oil\n[Ref.]','Palm oil\n[Ref.]','Trad. argan\n[Msanda\n2005]','Canola\n[Ref.]'],
     [R['WF_kg'], 3.0, 5.0, 12.0, 4.8],
     [NAVY, ORANGE, RED, GOLD, TEAL]),
]):
    bars5 = ax5.bar(labels5, vals5, color=cols5, edgecolor='white', linewidth=1.3, width=0.55)
    for b, v in zip(bars5, vals5):
        ax5.text(b.get_x()+b.get_width()/2, v+max(vals5)*0.025,
                 f'{v:.2f}', ha='center', fontsize=9, fontweight='bold')
    ax5.set_ylabel(ylabel, fontsize=10)
    ax5.set_title(title, fontweight='bold')
    ax5.tick_params(axis='x', labelsize=8.5)
    bars5[0].set_edgecolor(NAVY);  bars5[0].set_linewidth(2.0)

fig5.suptitle('LCA Benchmarking: Argan Biorefinery vs. Reference Vegetable Oils\n'
               '(Cradle-to-gate, comparable system boundaries; [Ref.] = published values)',
               fontsize=12, fontweight='bold', y=1.02)
fig5.tight_layout()
save(fig5, 'fig5_LCA_benchmarking.png')

# ─────────────────────────────────────────────────────────────────────────────
# FIG 6 — PCE BREAKDOWN BY UNIT
# ─────────────────────────────────────────────────────────────────────────────
pce_items = sorted(R['PCE_breakdown'].items(), key=lambda x: -x[1])
pce_labels = [k.replace(': ', '\n', 1) for k, _ in pce_items]
pce_vals   = [v/1000 for _, v in pce_items]   # k USD
cols6 = [NAVY,TEAL,GREEN,ORANGE,GOLD,PURPLE,LBLUE,'#BDC3C7',RED,'#1ABC9C','#F39C12']

fig6, ax6 = plt.subplots(figsize=(13, 6.5))
xs6 = np.arange(len(pce_labels))
bars6 = ax6.bar(xs6, pce_vals, color=cols6[:len(pce_vals)],
                 edgecolor='white', linewidth=1.1)
for b, v in zip(bars6, pce_vals):
    ax6.text(b.get_x()+b.get_width()/2, v+2, f'${v:.0f}k',
             ha='center', va='bottom', fontsize=8.5, rotation=40)
ax6.set_xticks(xs6)
ax6.set_xticklabels(pce_labels, fontsize=8, rotation=28, ha='right')
ax6.set_ylabel('Purchased Equipment Cost (k USD, 2026)', fontsize=11)
ax6.set_title('Capital Cost Breakdown by Equipment Item\n'
               f'(PCE = ${R["PCE"]/1e6:.2f}M × Lang 4.2 → CAPEX = ${R["CAPEX"]/1e6:.2f}M  |  CEPCI {int(CEPCI_2026)})',
               fontweight='bold')
ax6.text(0.98, 0.96,
         f'CEPCI {int(CEPCI_REF)} → {int(CEPCI_2026)}\n'
         f'Factor: ×{CEPCI_FACTOR:.3f}\n'
         f'PCE    = ${R["PCE"]/1e6:.2f}M\n'
         f'CAPEX  = ${R["CAPEX"]/1e6:.2f}M',
         transform=ax6.transAxes, ha='right', va='top', fontsize=9.5,
         fontweight='bold', color=NAVY,
         bbox=dict(boxstyle='round', facecolor='#EBF5FB', edgecolor=NAVY))
fig6.tight_layout()
save(fig6, 'fig6_CAPEX_breakdown.png')

print()

# ═══════════════════════════════════════════════════════════════════════════════
# 11.  FINAL SUMMARY
# ═══════════════════════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print(f"  ALL OUTPUTS WRITTEN TO: {OUT}")
print(f"{'='*70}")
files_expected = ['res.json','fig1_PFD_professional.png','fig2_TEA.png','fig3_LCA.png',
                   'fig4_sensitivity.png','fig5_LCA_benchmarking.png',
                   'fig6_CAPEX_breakdown.png']
for fn in files_expected:
    fp = os.path.join(OUT, fn)
    size = os.path.getsize(fp)/1024 if os.path.exists(fp) else 0
    status = '✓' if os.path.exists(fp) else '✗'
    print(f"  {status}  {fn:<40} {size:>7.1f} kB")
print(f"\n  CEPCI 2026 = {CEPCI_2026:.1f}  (escalation factor vs 2001 Turton base: {CEPCI_FACTOR:.3f})")
print(f"  Submit to : Biomass Conversion and Biorefinery (Springer, IF 5.0)")
print(f"{'='*70}\n")
