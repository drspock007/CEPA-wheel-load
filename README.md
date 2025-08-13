# CEPA Surface Loading – Mini App (Streamlit)

Une petite app web pour calculer les contraintes sur conduites enterrées soumises à une charge de surface.
MVP: charge ponctuelle au-dessus de la génératrice supérieure (superposition des empreintes roues/chenilles à venir).

## Lancer en ligne (le plus simple)
1. Crée un compte gratuit sur **Streamlit Community Cloud**.
2. Crée un nouveau projet en important ces fichiers (`app.py`, `cepa_surface_loading.py`, `requirements.txt`).
3. Set `app.py` comme fichier principal. Déploie — c’est tout.

## Lancer en local (si tu veux tester sur ton Mac)
```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

## Fichiers
- `app.py` : l'app Streamlit
- `cepa_surface_loading.py` : moteur de calcul
- `requirements.txt` : dépendances Python

## Note
- Paramètres sol CEPA (Kb, Kz) exposés pour calibration.
- Prochaine itération : empreintes véhicules (roues/chenilles), pass/fail codes (B31.4/B31.8/CSA), & rapport PDF.
