"""Generate lib/tests/wrf_templates.esm: every WRF template applied to fixed
inputs, with expected values from independent binary64 Python formulas
written from the Fortran (not from the JSON bodies).

Usage: python3 lib/tests/gen_wrf_templates.py lib/tests/wrf_templates.esm
"""
import json, math, sys

# WRF constants (lib/wrf_constants.esm)
r_d, r_v, cliq, cice, psat = 287.0, 461.6, 4190.0, 2106.0, 610.78
p1000mb, xlv, xls, svpt0 = 100000.0, 2.5e6, 2.85e6, 273.15
svp1, svp2, svp3 = 0.6112, 17.67, 29.65
cp = 7 * r_d / 2; cv = cp - r_d; cpv = 4 * r_v
ep_1 = r_v / r_d - 1; ep_2 = r_d / r_v; rcp = r_d / cp
cpovcv = cp / cv; rvovrd = r_v / r_d; svp1_pa = 1000 * svp1

# inputs
IN = dict(T=283.0, p=85000.0, qv=0.008, qc=2.0e-4, qi=5.0e-5, theta=295.0,
          ttp=273.16, theta_m=303.0, alpha_d=0.95, rho=1.1)
T, p, qv, qc, qi, theta = (IN[k] for k in ("T", "p", "qv", "qc", "qi", "theta"))
ttp, theta_m, alpha_d, rho = (IN[k] for k in ("ttp", "theta_m", "alpha_d", "rho"))

V = {}   # name -> (expected, units, description, rhs)
def tpl(name, **b):
    return {"op": "apply_expression_template", "args": [], "name": name, "bindings": b}
def c(value, units):
    """A literal that declares its units (esm 1.2.0, spec 4.8.5)."""
    return {"op": "const", "args": [], "value": value, "units": units}
def add(var, expected, units, desc, rhs):
    V[var] = (expected, units, desc, rhs)

exner = (p / p1000mb) ** rcp                     # pi = (p/p1000mb)**rcp
add("exner", exner, "1", "exner_function at p with p0 = p1000mb, kappa = rcp.",
    tpl("exner_function", p="p", p0="wrf.p1000mb", kappa="wrf.rcp"))
add("T_from_theta", theta * exner, "K", "temperature_from_theta: theta pi.",
    tpl("temperature_from_theta", theta="theta", exner="exner"))
add("theta_roundtrip", theta, "K", "theta_from_temperature(temperature_from_theta(theta)) recovers theta.",
    tpl("theta_from_temperature", T="T_from_theta", exner="exner"))
add("theta_of_T", T / exner, "K", "theta_from_temperature: T / pi.",
    tpl("theta_from_temperature", T="T", exner="exner"))
add("theta_v", theta * (1 + ep_1 * qv), "K", "virtual_potential_temperature.",
    tpl("virtual_potential_temperature", theta="theta", qv="qv", ep_1="wrf.ep_1"))
# YSU: thlix = (tx - xlv*qcx/cp - xls*qix/cp)/pi2d with tx = thx*pi2d
add("theta_li", (theta * exner - xlv * qc / cp - xls * qi / cp) / exner, "K",
    "liquid_ice_potential_temperature, against YSU's (tx - xlv qc/cp - xls qi/cp)/pi.",
    tpl("liquid_ice_potential_temperature", theta="theta", qc="qc", qi="qi", exner="exner",
        xlv="wrf.xlv", xls="wrf.xls", cp="wrf.cp"))
add("rho_moist", p / (r_d * T * (1 + ep_1 * qv)), "kg/m^3", "moist_air_density.",
    tpl("moist_air_density", p="p", T="T", qv="qv", r_d="wrf.r_d", ep_1="wrf.ep_1"))
# WSM6 fpvs, water branch
xa = -(cpv - cliq) / r_v
xb = xa + xlv / (r_v * ttp)
add("xa", xa, "1", "fpvs_exponent_a over water: -(cpv - cliq)/r_v.",
    tpl("fpvs_exponent_a", c_vapor="wrf.cpv", c_condensate="wrf.cliq", r_v="wrf.r_v"))
add("xb", xb, "1", "fpvs_exponent_b over water: xa + xlv/(r_v ttp).",
    tpl("fpvs_exponent_b", xa="xa", L="wrf.xlv", r_v="wrf.r_v", ttp="ttp"))
tr = ttp / T
es_fpvs = psat * math.exp(math.log(tr) * xa) * math.exp(xb * (1 - tr))  # WSM6: psat*exp(log(tr)*xa)*exp(xb*(1-tr))
add("es_fpvs", es_fpvs, "Pa", "saturation_vapor_pressure_fpvs over water at T, against WSM6's psat exp(log(tr) xa) exp(xb (1 - tr)), tr = ttp/T.",
    tpl("saturation_vapor_pressure_fpvs", T="T", psat="wrf.psat", ttp="ttp", xa="xa", xb="xb"))
add("es_fpvs_at_ttp", psat, "Pa", "saturation_vapor_pressure_fpvs at T = ttp equals psat, the anchor of the fit.",
    tpl("saturation_vapor_pressure_fpvs", T="ttp", psat="wrf.psat", ttp="ttp", xa="xa", xb="xb"))
add("qsat_fpvs", ep_2 * es_fpvs / (p - es_fpvs), "kg/kg", "saturation_mixing_ratio from es_fpvs.",
    tpl("saturation_mixing_ratio", es="es_fpvs", p="p", ep_2="wrf.ep_2"))
xlv1 = cliq - cpv                                # WSM6 xlv1 = cl - cpv; xlcal = xlv0 - xlv1 (t - t0c), xlv0 = xlv
add("xl", xlv - xlv1 * (T - svpt0), "J/kg", "latent_heat_of_vaporization with WSM6's bindings (T_ref = svpt0).",
    tpl("latent_heat_of_vaporization", T="T", xlv="wrf.xlv", cliq="wrf.cliq", cpv="wrf.cpv", T_ref="wrf.svpt0"))
add("xl_at_ref", xlv, "J/kg", "latent_heat_of_vaporization at T = T_ref equals xlv.",
    tpl("latent_heat_of_vaporization", T="wrf.svpt0", xlv="wrf.xlv", cliq="wrf.cliq", cpv="wrf.cpv", T_ref="wrf.svpt0"))
add("cpm", cp * (1 - qv) + qv * cpv, "J/(kg*K)", "moist_heat_capacity (WSM6 cpmcal).",
    tpl("moist_heat_capacity", qv="qv", cp="wrf.cp", cpv="wrf.cpv"))
es_b = svp1_pa * math.exp(svp2 * (T - svpt0) / (T - svp3))
add("es_bolton", es_b, "Pa", "bolton_saturation_vapor_pressure with svp1_pa, so in Pa.",
    tpl("bolton_saturation_vapor_pressure", T="T", svp1="wrf.svp1_pa", svp2="wrf.svp2", svp3="wrf.svp3", svpt0="wrf.svpt0"))
add("qsat_bolton", ep_2 * es_b / (p - es_b), "kg/kg", "saturation_mixing_ratio from es_bolton with p in Pa (sfclayrev's form, with pressures in Pa rather than kPa).",
    tpl("saturation_mixing_ratio", es="es_bolton", p="p", ep_2="wrf.ep_2"))
# YSU rvls = 100.*6.112*exp(17.67*(t-273.16)/(t-29.65))*(ep2/p)
add("rvls", 100. * 6.112 * math.exp(17.67 * (T - 273.16) / (T - 29.65)) * (ep_2 / p), "kg/kg",
    "bolton_saturation_mixing_ratio against bl_ysu.F90's rvls with T0 = 273.16 and es0 = 611.2 Pa.",
    tpl("bolton_saturation_mixing_ratio", T="T", T0=c(273.16, "K"), p="p", es0=c(611.2, "Pa"), a=c(17.67, "1"), b=c(29.65, "K"), ep_2="wrf.ep_2"))
add("rho_dry", p / (r_d * T), "kg/m^3", "dry_air_density.",
    tpl("dry_air_density", p="p", T="T", r_d="wrf.r_d"))
add("theta_dry", theta_m / (1 + rvovrd * qv), "K", "moist_theta_to_dry_theta.",
    tpl("moist_theta_to_dry_theta", theta_m="theta_m", qv="qv", rvovrd="wrf.rvovrd"))
# calc_p_rho_phi: p = p1000mb*((r_d*theta_m)/(p1000mb*alpha_d))**cpovcv
add("p_eos", p1000mb * ((r_d * theta_m) / (p1000mb * alpha_d)) ** cpovcv, "Pa", "wrf_equation_of_state against calc_p_rho_phi.",
    tpl("wrf_equation_of_state", theta_m="theta_m", alpha_d="alpha_d", p0="wrf.p1000mb", r_d="wrf.r_d", cpovcv="wrf.cpovcv"))
add("rho_total", (1 + qv) / alpha_d, "kg/m^3", "moist_density_from_inverse_dry_density.",
    tpl("moist_density_from_inverse_dry_density", qv="qv", alpha_d="alpha_d"))
# WSM6 statement functions
nu = 1.496e-6 * T ** 1.5 / (T + 120.) / rho
add("nu", nu, "m^2/s", "kinematic_viscosity_of_air against WSM6 viscos(t, den) = 1.496e-6 t**1.5/(t + 120)/den.",
    tpl("kinematic_viscosity_of_air", T="T", rho="rho", coef=c(1.496e-6, "kg/(m*s*K^(1/2))"), T_sutherland=c(120.0, "K")))
add("Dv", 8.794e-5 * T ** 1.81 / p, "m^2/s", "diffusivity_of_water_vapor against WSM6 diffus(t, p) = 8.794e-5 t**1.81/p.",
    tpl("diffusivity_of_water_vapor", T="T", p="p", coef=c(8.794e-5, "m^2*Pa/(s*K^1.81)"), expn=1.81))
add("ka", 1.414e3 * nu * rho, "J/(m*s*K)", "thermal_conductivity_of_air against WSM6 xka(t, den) = 1.414e3 viscos(t, den) den.",
    tpl("thermal_conductivity_of_air", nu="nu", rho="rho", coef=c(1.414e3, "J/(kg*K)")))

units_in = dict(T="K", p="Pa", qv="kg/kg", qc="kg/kg", qi="kg/kg", theta="K", ttp="K",
                theta_m="K", alpha_d="m^3/kg", rho="kg/m^3")
variables = {k: {"type": "parameter", "units": units_in[k], "default": v,
                 "description": f"Test input {k}."} for k, v in IN.items()}
for k, (e, u, d, _) in V.items():
    variables[k] = {"type": "unknown", "units": u, "description": d}
doc = {
  "esm": "1.2.0",
  "metadata": {
    "name": "WRFTemplateChecks",
    "description": "Checks for the expression templates of lib/wrf_thermo.esm and lib/wrf_air_properties.esm. A template-library file cannot carry inline tests, so this 0-D component mounts lib/wrf_constants.esm, imports both libraries and applies every template once to fixed inputs (T = 283 K, p = 850 hPa, qv = 8 g/kg, ...), binding the constants the way the WRF schemes do. The expected values are the WRF Fortran expressions (WSM6 fpvs, xlcal, cpmcal, viscos, diffus, xka; YSU thlix and rvls; calc_p_rho_phi; Bolton e_s) written out independently in Python and evaluated in binary64, plus the identities that anchor each fit: fpvs(ttp) = psat, L(T_ref) = xlv and theta -> T -> theta. Every literal bound to a template carries its units (esm 1.2.0 const), so each equation is dimensionally checked. Tolerance 1e-10 relative, to allow libm differences in exp, log and pow.",
    "tags": ["wrf", "thermodynamics", "expression-templates", "tests"],
    "created": "2026-10-15T00:00:00Z",
    "authors": ["EqWeFiC"]
  },
  "models": {
    "WRFTemplateChecks": {
      "expression_template_imports": [
        {"ref": "../wrf_thermo.esm"},
        {"ref": "../wrf_air_properties.esm"}
      ],
      "subsystems": {"wrf": {"ref": "../wrf_constants.esm"}},
      "variables": variables,
      "equations": [{"lhs": k, "rhs": r} for k, (_, _, _, r) in V.items()],
      "tests": [{
        "id": "templates_match_fortran",
        "description": "Every template against its WRF Fortran expression evaluated in binary64 (generator: lib/tests/gen_wrf_templates.py).",
        "time_span": {"start": 0.0, "end": 1.0},
        "tolerance": {"rel": 1e-10},
        "assertions": [{"variable": k, "time": 0.0, "expected": e} for k, (e, _, _, _) in V.items()]
      }]
    }
  }
}
open(sys.argv[1], "w").write(json.dumps(doc, indent=2, ensure_ascii=False) + "\n")
