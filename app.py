
import streamlit as st
import pandas as pd
from math import isnan
from dataclasses import asdict
from typing import Optional
import os

# Import engine
import sys
sys.path.append(os.path.dirname(__file__))
from cepa_surface_loading import (
    Pipe, Soil, Misc, example_compute_over_crown
)

st.set_page_config(page_title="CEPA Surface Loading (MVP)", page_icon="🛠️", layout="wide")

st.title("CEPA Surface Loading – Mini App (MVP)")
st.caption("Prototype simplifié : charge ponctuelle sur la génératrice supérieure (superposition véhicules à venir).")

with st.sidebar:
    st.header("Unités & méthode")
    units = st.radio("Système d'unités", ["SI", "IMP"], index=0, help="SI (m, Pa) ou Impérial (in, ft, psi).")
    eqv = st.radio("Contrainte équivalente", ["tresca","von_mises"], index=0)
    vehicle_type = st.selectbox("Type de véhicule", ["highway","farm_construction","track"])
    pavement = st.selectbox("Type de revêtement", ["flexible","rigid"])
    kb = st.number_input("Kb (paramètre sol)", value=0.031, step=0.001, format="%.3f")
    kz = st.number_input("Kz (paramètre sol)", value=0.0915, step=0.001, format="%.4f")
    st.markdown("> **Astuce**: Garde Tresca et un revêtement *flexible* pour rester conservateur.")

st.subheader("1) Conduite")
col1, col2, col3 = st.columns(3)
with col1:
    D = st.number_input("Diamètre extérieur D", value=0.3239 if units=="SI" else 12.75, min_value=0.0)
with col2:
    t = st.number_input("Épaisseur t", value=0.00953 if units=="SI" else 0.375, min_value=0.0)
with col3:
    SMYS = st.number_input("SMYS", value=359e6 if units=="SI" else 52000.0, min_value=0.0)

col1, col2, col3 = st.columns(3)
with col1:
    MOP = st.number_input("MOP (pression d'exploitation max.)", value=7200000.0 if units=="SI" else 1044.0, min_value=0.0)
with col2:
    dT = st.number_input("ΔT (Tinstallé - Topération)", value=11.1 if units=="SI" else 20.0)
with col3:
    theta = st.number_input("Angle de lit θ (deg)", value=30.0, min_value=0.0, max_value=180.0)

st.markdown("---")
st.subheader("2) Sol")
col1, col2, col3 = st.columns(3)
with col1:
    rho = st.number_input("Densité ρ", value=1920.0 if units=="SI" else 120.0, min_value=0.0)
with col2:
    H = st.number_input("Profondeur de recouvrement H", value=2.0 if units=="SI" else 6.0, min_value=0.0)
with col3:
    Eprime = st.number_input("Module de réaction E′", value=69000.0 if units=="SI" else 10.0, min_value=0.0,
                              help="SI: Pa ; Impérial: psi")

st.markdown("---")
st.subheader("3) Charge ponctuelle (MVP)")
col1, col2 = st.columns(2)
with col1:
    F_point = st.number_input("Charge ponctuelle au sol (par point)", value=20000.0 if units=="SI" else 4500.0, min_value=0.0,
                              help="Approche conservative: répartir la charge roue/chenille en points et prendre le plus défavorable (sur la génératrice).")
with col2:
    st.write("**Impact factor** selon véhicule & revêtement appliqué automatiquement (peut être raffiné dans une prochaine version).")

# Build dataclasses
pipe = Pipe(D=D, t=t, MOP=MOP, SMYS=SMYS, dT=dT, units=units)
soil = Soil(rho=rho, H=H, theta_deg=theta, Eprime=Eprime, units=units)
misc = Misc(vehicle_type=vehicle_type, pavement=pavement, eqv=eqv)

if st.button("Calculer"):
    try:
        res = example_compute_over_crown(pipe=pipe, soil=soil, misc=misc, F_point_surface=F_point, Kb=kb, Kz=kz)
        # Display results
        st.success("Calcul terminé ✅")
        # Pretty table
        df = pd.DataFrame({
            "Grandeur":[
                "P_live (Pa)","P_soil (Pa)",
                "σH_int (Pa)","σH_live_bend (Pa)","σH_soil_bend (Pa)","σH_total (Pa)",
                "σL_beam (Pa)","σL_total (Pa)",
                "σ_eq (Pa)","SMYS (Pa)","Utilisation σ_eq/SMYS"
            ],
            "Valeur":[
                res["P_live_Pa"],res["P_soil_Pa"],
                res["sigma_H_int_Pa"],res["sigma_H_live_bend_Pa"],res["sigma_H_soil_bend_Pa"],res["sigma_H_total_Pa"],
                res["sigma_L_beam_Pa"],res["sigma_L_total_Pa"],
                res["sigma_eq_Pa"],res["SMYS_Pa"],res["utilization_eq"]
            ]
        })
        st.dataframe(df, use_container_width=True)

        # Download Excel (inputs + outputs)
        xls_name = "resultats_cepa_mvp.xlsx"
        with pd.ExcelWriter(xls_name) as writer:
            # inputs
            pd.DataFrame([pipe.__dict__]).to_excel(writer, index=False, sheet_name="Pipe")
            pd.DataFrame([soil.__dict__]).to_excel(writer, index=False, sheet_name="Soil")
            pd.DataFrame([misc.__dict__]).to_excel(writer, index=False, sheet_name="Misc")
            df.to_excel(writer, index=False, sheet_name="Results")
        with open(xls_name, "rb") as f:
            st.download_button("📥 Télécharger Excel (entrées + résultats)", f, file_name=xls_name)

        st.info("Note: Ce MVP calcule pour 1 point de charge sur la génératrice. Les empreintes complètes (roues/chenilles/grille) arrivent à la prochaine itération.")
    except Exception as e:
        st.error(f"Erreur de calcul: {e}")
