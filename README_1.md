# Capilaridad y Flujo en Xilema
### Simulación computacional del transporte de agua en plantas vasculares

**Autores:** Lucas Neisskenwirth · Martín Costa · Max Lavandero  
**Curso:** Universidad Adolfo Ibáñez, Santiago, 2026

---

## Resumen

Este proyecto simula dos fenómenos físicos que explican el transporte de agua en el xilema de las plantas:

1. **Capilaridad (Ley de Jurin):** ascenso del agua en tubos estrechos en función del radio.
2. **Flujo viscoso (Hagen-Poiseuille):** caudal en red de vasos bajo gradiente de potencial hídrico.

El modelo de red de poros resuelve el sistema de ecuaciones de Kirchhoff para obtener perfiles de potencial hídrico y caudal en estado estacionario.

---

## Arquitectura del modelo

```
Raíz (ψ = 0)
    │
    ├─[seg 1, r₁, K₁]─ nodo 1
    ├─[seg 2, r₂, K₂]─ nodo 2
    │    ...
    └─[seg N, rₙ, Kₙ]─ Hoja (ψ = ΔΨ)

Kᵢ = π rᵢ⁴ / (8 η L)   [Hagen-Poiseuille]
ΣQ en nodo interno = 0  [Kirchhoff]
→ Sistema lineal A·ψ = b resuelto con scipy.linalg.solve
```

---

## Estructura del repositorio

```
xilema_simulacion/
├── README.md
├── requirements.txt
├── codigo/
│   ├── simulacion_xilema.py   ← script principal
│   └── generar_informe.js     ← genera el informe Word
└── outputs/
    ├── fig1_jurin_ascenso.png
    ├── fig2_hagen_poiseuille.png
    ├── fig3_red_poros_base.png
    ├── fig4_sensibilidad_heatmap.png
    ├── fig5_comparacion_configs.png
    ├── tabla_resultados.csv
    └── Informe_Xilema_UAI.docx
```

---

## Requisitos e instalación

**Python 3.9+**

```bash
pip install -r requirements.txt
```

`requirements.txt` incluye: `numpy`, `scipy`, `matplotlib`

---

## Cómo ejecutar

```bash
# 1. Clonar / descargar el repo
cd xilema_simulacion/codigo

# 2. Ejecutar la simulación (genera todos los gráficos y la tabla)
python simulacion_xilema.py

# Los archivos de salida se guardan en ../outputs/
```

---

## Experimentos generados

| Figura | Descripción |
|--------|-------------|
| fig1 | Ley de Jurin: altura capilar vs radio (1–200 µm) |
| fig2 | Hagen-Poiseuille: caudal vs radio (ΔP = 1 MPa) |
| fig3 | Red de poros base: perfil de ψ y caudal en estado estacionario |
| fig4 | Mapa de calor sensibilidad: caudal vs (r, ΔΨ) |
| fig5 | Comparación anatómica: 4 configuraciones de vasos |

---

## Resultados esperados

Al ejecutar la simulación obtendrás:

- **Jurin (r = 20 µm):** h ≈ 74 cm — insuficiente para árboles altos.
- **Hagen-Poiseuille:** duplicar r multiplica Q por 16.
- **Red base:** perfil lineal de ψ, caudal constante en estado estacionario.
- **Sensibilidad:** el radio domina el caudal con sensibilidad r⁴; ΔΨ actúa linealmente.

---

## Semilla aleatoria

```python
np.random.seed(42)
```
Los resultados son **completamente deterministas** (no hay componente estocástico en el modelo base).

---

## Referencias

- Koch et al. (2004). *Nature*, 428, 851–854.
- Allen et al. (2010). *Forest Ecology and Management*, 259, 660–684.
- Tyree & Zimmermann (2002). *Xylem Structure and the Ascent of Sap*. Springer.
- Dixon & Joly (1895). *Phil. Trans. Royal Soc. B*, 186, 563–576.

---

## Licencia

Proyecto académico UAI 2026. Uso libre con atribución.
