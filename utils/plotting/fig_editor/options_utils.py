import matplotlib as mpl

NOLEGEND = '_nolegend_'

def get_linestyles():
    return {'-': 'Solid',
            '--': 'Dashed',
            '-.': 'DashDot',
            ':': 'Dotted',
            'None': 'None',
            }


def get_drawstyles():
    return {
        'default': 'Default',
        'steps-pre': 'Steps (Pre)',
        'steps': 'Steps',
        'steps-mid': 'Steps (Mid)',
        'steps-post': 'Steps (Post)'}


def get_ticks_directions():
    return ["in", "out", "inout"]


def get_scales():
    return ['linear', 'log', 'symlog', 'logit']


def regenerate_legend(axes):
    draggable = None
    loc = None
    ncol = 1
    old_legend = axes.get_legend()
    if old_legend is not None:
        draggable = old_legend._draggable is not None
        ncol = old_legend._ncol
        loc = old_legend._get_loc()
    new_legend = axes.legend(ncol=ncol)
    if new_legend:
        new_legend.draggable(draggable)
        if loc:
            new_legend._set_loc(loc)
    return new_legend


def get_markers():
    return mpl.markers.MarkerStyle.markers


def is_errorbar(o):
    return isinstance(o, mpl.container.ErrorbarContainer)


def apply_to_ebar(func, ebar,  value):
    """ Do function `func` with input `line, value` to all lines of the errorbar ebar """
    func(ebar[0], value)
    for collection in ebar[1:]:
        for bar in collection:
            func(bar, value)


def get_errorbar_from_line(line):
    """ Returns first container of type errorbar if line is in there, otherwise index errror. """
    return [cont for cont in line.axes.containers if is_errorbar(cont) and line in cont][0]


def prepare_formdata(d, init):
    """Prepare entry for FormLayout.

    `d` is a mapping of shorthands to style names (a single style may
    have multiple shorthands, in particular the shorthands `None`,
    `"None"`, `"none"` and `""` are synonyms); `init` is one shorthand
    of the initial style.

    This function returns a tuple of the initial value and the choices of the formdata:
    `initial_name, [(shorthand, style_name), (shorthand, style_name), ...]`.
    """
    # Drop duplicate shorthands from dict (by overwriting them during
    # the dict comprehension).
    name2short = {name: short for short, name in d.items()}
    # Convert back to {shorthand: name}.
    short2name = {short: name for name, short in name2short.items()}
    # Find the kept shorthand for the style specified by init.
    canonical_init = name2short[d[init]]
    # Sort by representation and prepend the initial value.
    return canonical_init, sorted(short2name.items(), key=lambda short_and_name: short_and_name[1])
