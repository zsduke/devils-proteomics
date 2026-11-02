"""Shared publication graphics and self-contained, offline mathematical typesetting."""
from copy import deepcopy
from functools import lru_cache
from pathlib import Path
import base64
import html
import re

COLORS = {"up": "#D55E00", "down": "#0072B2", "neutral": "#B3B3B3", "unchanged": "#555555"}
CATEGORICAL = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000"]
DIVERGING = [[0, "#0072B2"], [0.5, "#FFFFFF"], [1, "#D55E00"]]
FIGURE_REFERENCES = [
    {"title": "Nature research figure guide", "url": "https://research-figure-guide.nature.com/figures/preparing-figures-our-specifications/"},
    {"title": "Ctortecka et al., Nature Communications 15, 5707 (2024)", "url": "https://doi.org/10.1038/s41467-024-49651-w"},
    {"title": "Wong, Colour blindness, Nature Methods 8, 441 (2011)", "url": "https://doi.org/10.1038/nmeth.1618"},
]


def publication_layout(title="", **kwargs):
    axis = {"showline": True, "linecolor": "#111111", "linewidth": 1,
            "ticks": "outside", "ticklen": 5, "tickwidth": 1, "tickcolor": "#111111",
            "showgrid": False, "zeroline": False, "mirror": False, "automargin": True}
    layout = {"template": "none", "title": {"text": title, "x": 0.02, "xanchor": "left"},
              "paper_bgcolor": "#FFFFFF", "plot_bgcolor": "#FFFFFF",
              "font": {"family": "Arial", "size": 13, "color": "#111111"},
              "colorway": list(CATEGORICAL), "margin": {"t": 65, "b": 75, "l": 75, "r": 25},
              "legend": {"borderwidth": 0, "font": {"color": "#111111"}, "bgcolor": "rgba(255,255,255,0)"}}
    for key in ("xaxis", "yaxis"):
        layout[key] = {**axis, **kwargs.pop(key, {})}
    layout.update(deepcopy(kwargs))
    return layout


def math_markup(tex, display=False):
    """Escape markup; callers supply formula source, never executable HTML."""
    tag = "div" if display else "span"
    return f'<{tag} class="math-formula" data-display="{str(display).lower()}" data-tex="{html.escape(tex, quote=True)}">{html.escape(tex)}</{tag}>'


@lru_cache(maxsize=1)
def math_assets():
    """Embed the pinned KaTeX runtime and WOFF2 fonts so exported HTML is portable."""
    vendor = Path(__file__).parent / "static/vendor/katex"
    css = (vendor / "katex.min.css").read_text("utf-8")
    def inline_font(match):
        filename = re.search(r'url\((fonts/[^)]+\.woff2)\)', match.group())
        if not filename:
            raise ValueError("Missing WOFF2 font in bundled KaTeX stylesheet.")
        data = base64.b64encode((vendor / filename.group(1)).read_bytes()).decode("ascii")
        return 'src:url(data:font/woff2;base64,' + data + ') format("woff2")'
    css = re.sub(r'src:[^;}]+', inline_font, css)
    script = (vendor / "katex.min.js").read_text("utf-8").replace("</script", "<\\/script")
    notices = "<!--\nEmbedded third-party licenses\n" + (vendor / "LICENSE").read_text("utf-8") + "\n" + (Path(__file__).parent / "licenses/PLOTLY_MIT.txt").read_text("utf-8") + "\n-->"
    return (notices + '<style>' + css + '\n.math-formula[data-display="true"]{overflow-x:auto;padding:6px 0}'
            '.math-formula{color:#111}.katex-display{margin:.65em 0}</style><script>' + script + '</script>'
            '<script>document.querySelectorAll("[data-tex]").forEach(function(el){'
            'katex.render(el.dataset.tex,el,{displayMode:el.dataset.display==="true",'
            'throwOnError:false,trust:false,strict:"error",output:"htmlAndMathml"});});</script>')


def latex_number(value):
    text = format(value, ".12g")
    if "e" in text:
        mantissa, exponent = text.split("e")
        return mantissa + r"\times 10^{" + str(int(exponent)) + "}"
    return text


def volcano_formulae(config):
    bx, by = latex_number(config["x_base"]), latex_number(config["y_base"])
    metric = config["significance_metric"]
    equations = [rf'x=\log_{{{bx}}}(\mathrm{{FC}})',
                 rf'y=-\log_{{{by}}}\!\left(\max\!\left({metric},10^{{-300}}\right)\right)']
    cutoff = latex_number(config["fc_cutoff"])
    if config["fc_cutoff_mode"] == "log":
        effect = rf'|x|\geq {cutoff}'
    else:
        effect = rf'\mathrm{{FC}}\geq {cutoff}\quad\mathrm{{or}}\quad\mathrm{{FC}}\leq \frac{{1}}{{{cutoff}}}'
    alpha = latex_number(config["significance_cutoff"])
    significance = (rf'{metric}\leq {alpha}' if config["significance_mode"] == "probability"
                    else rf'-\log_{{{by}}}({metric})\geq {alpha}')
    equations.append(effect + r'\qquad\mathrm{and}\qquad ' + significance)
    return ''.join(math_markup(eq, display=True) for eq in equations)
