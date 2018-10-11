"""
TODO:
 Lines:
 - Everything

 Axes:
 - Everything

 - Add visibility to options_figure

"""
import matplotlib as mpl
import matplotlib.backends.qt_editor.formlayout as formlayout

import six

from gui_utils import get_icon
import options_utils as outils
from general_helper import suppress_exception


def change_properties(artist, parent=None):
    """ Change the properties of artist via user interface. """
    regen_legend = False
    if isinstance(artist, mpl.lines.Line2D):
        try:
            ebar = outils.get_errorbar_from_line(artist)
        except IndexError:
            props, title = _get_line_properties(artist)
        else:
            props, title, = _get_ebar_properties(ebar)
        regen_legend = True

    elif isinstance(artist, mpl.axis.Axis):
        props, title = _get_axis_properties(artist)

    elif isinstance(artist, mpl.text.Text):
        props, title = _get_text_properties(artist)

    elif isinstance(artist, mpl.legend.Legend):
        props, title = _get_legend_properties(artist)
        regen_legend = True

    # elif isinstance(artist, mpl.axes.Axes):
    #     props, title = _get_axes_properties(artist)

    else:
        return

    data = _show_form(props, title, parent)

    if data is not None:
        _apply_data_to_properties(props, data, parent)
        if regen_legend:
            new_leg = outils.regenerate_legend(artist.axes)


def _show_form(props, title, parent):
    form = []
    for ppty in props:
        if isinstance(ppty, Property):
            # proper property
            form.append(ppty.get_form_line())
        elif isinstance(ppty, six.string_types):
            # tuple
            form.append((None, "<b>{}</b>".format(ppty)))
        else:
            # separator
            form.append((None, None))

    return formlayout.fedit(form, title=title, parent=parent, icon=get_icon('edit'),
                            apply=
                            lambda d, p=props, par=parent: _apply_data_to_properties(p, d, par)
                            )


# Form layouts #################################################################


def _get_axes_properties(ax):
    title = "Edit Axes"
    props = ["Axes Edit Content coming up."]
    return props, title


def _get_axis_properties(axis):
    """ Properties for XAxis, YAxis """
    title = "Edit {}Axis".format(axis.axis_name.upper())

    plane, pos_choices = _get_plane_from_axis(axis)

    # defaults for spines
    spines = [axis.axes.spines[pos] for pos in pos_choices]
    color = spines[0].get_edgecolor()
    width = spines[0].get_linewidth()
    try:
        ticks_dir = axis.get_ticks_direction()[0]
    except IndexError:
        ticks_dir = outils.get_ticks_directions()[0]

    # defaults for grid
    grid_lines = axis.axes.__getattribute__("get_{}gridlines".format(plane))()
    if grid_lines:
        grid_width = grid_lines[0].get_linewidth()
        grid_color = grid_lines[0].get_color()
        grid_style = grid_lines[0].get_linestyle()
    else:
        grid_width = .8
        grid_color = '#b0b0b0'
        grid_style = outils.get_linestyles().keys()[0]

    grid_style_data = outils.prepare_formdata(outils.get_linestyles(), grid_style)

    # find limits
    lim = axis.axes.__getattribute__("get_{}lim".format(plane))()
    set_lim_func = axis.axes.__getattribute__("set_{}lim".format(plane))

    props = [
        "Labels:",
        Property("Label", axis.set_label_text,
                 axis.get_label_text(),
                 None),
        Property("Label Position",
                 axis.set_label_position, axis.get_label_position(),
                 pos_choices, ),
        Property("Label Padding",
                 lambda x: setattr(axis, "labelpad", x), axis.labelpad,
                 None),
        "Spine:",
        Property("Scale",
                 axis._set_scale, axis.get_scale(),
                 outils.get_scales()),
        Property("Lower Limit",
                 lambda x: set_lim_func(x, None), lim[0],
                 None),
        Property("Upper Limit",
                 lambda x: set_lim_func(None, x), lim[1],
                 None),
        Property("Tick Position",
                 lambda x, a=axis, p=plane:
                 axis.axes.tick_params(axis=p, which="both", direction=x),
                 ticks_dir,
                 outils.get_ticks_directions()),
        Property("Major Ticks",
                 lambda x, a=axis: _show_ticks(x, a, which="major"),
                 bool(axis.get_major_ticks()),
                 None),
        Property("Minor Ticks",
                 lambda x, a=axis: _show_ticks(x, a, which="minor"),
                 bool(axis.get_minor_ticks()),
                 None),
        Property("Color",
                 lambda x, a=axis: _set_axis_colors(x, a), color,
                 None),
        Property("Width",
                 lambda x: [s.set_linewidth(x) for s in spines], width,
                 None),
        "Grid:",
        Property("Line",
                 lambda x: axis.grid(ls=x), grid_style_data[0],
                 grid_style_data[1:]),
        Property("Width",
                 lambda x: axis.grid(lw=x), grid_width,
                 None),
        Property("Color",
                 lambda x: axis.grid(color=x), grid_color,
                 None),
        Property("Major",
                 lambda x: axis.grid(x, which='major'), axis._gridOnMajor,
                 None),
        Property("Minor",
                 lambda x: axis.grid(x, which='minor'), axis._gridOnMinor,
                 None),
        # Property("Visible",
        #          axis.grid, axis._gridOnMinor or axis._gridOnMajor,
        #          None),
    ]

    return props, title


def _get_line_properties(line):
    title = "Edit Line '{}'".format(line.get_label())
    ls_def, ls_choices = outils.prepare_formdata(outils.get_linestyles(), line.get_linestyle())
    ds_def, ds_choices = outils.prepare_formdata(outils.get_drawstyles(), line.get_drawstyle())
    ms_def, ms_choices = outils.prepare_formdata(outils.get_markers(), line.get_marker())

    props = [
        Property("Label",
                 line.set_label, line.get_label(),
                 None),
        Property("Visible",
                 line.set_visible, line.get_visible(),
                 None),
        Property("Z-Order",
                 line.set_zorder, line.get_zorder(),
                 None),
        "Line:",
        Property("Linestyle",
                 line.set_linestyle, ls_def,
                 ls_choices),
        Property("Drawstyle",
                 line.set_drawstyle, ds_def,
                 ds_choices),
        Property("Width",
                 line.set_linewidth, line.get_linewidth(),
                 None),
        Property("Color",
                 line.set_color, line.get_color(),
                 None),
        "Marker:",
        Property("Style",
                 line.set_marker, ms_def,
                 ms_choices),
        Property("Size",
                 line.set_markersize, line.get_markersize(),
                 None),
        Property("Edge Color",
                 line.set_markeredgecolor, line.get_markeredgecolor(),
                 None),
        Property("Edge Size",
                 line.set_markeredgewidth, line.get_markeredgewidth(),
                 None),
        Property("Face Color",
                 line.set_markerfacecolor, line.get_markerfacecolor(),
                 None),
    ]
    return props, title


def _get_ebar_properties(ebar):
    title = "Edit Line '{}'".format(ebar.get_label())

    bar = None  # bar reference
    with suppress_exception(IndexError):
        bar = ebar[2][0]

    label = ebar.get_label()
    if label is None:
        label = ""  # workaround as ebars can return None labels

    props = [
        Property("Label",
                 ebar.set_label, label,
                 None),
        Property("Visible",
                 lambda x, ebar=ebar: outils.apply_to_ebar(lambda e, v: e.set_visible(v), ebar, x),
                 ebar[0].get_visible(),
                 None),
        Property("Z-Order",
                 lambda x, ebar=ebar: outils.apply_to_ebar("set_zorder", ebar, x),
                 ebar[0].get_zorder(),
                 None),
    ]
    props += _get_line_properties(ebar[0])[0][3:]

    if bar is not None:
        ls = bar.get_linestyles()[0][0]
        if ls is None:
            ls = "None"  # it seems as if None here is "Solid" for some reason.
            # Workaround: set width to 0

        ls_def, ls_choices = outils.prepare_formdata(outils.get_linestyles(), ls)

        # ms_def, ms_choices = outils.prepare_formdata(
        #     outils.get_markers(), ebar[1][0].get_marker())
        props += [
            "Error Bars",
            Property("Linestyle",
                     lambda x, ebar=ebar: outils.apply_to_ebar(
                         "set_linestyle", ebar, x, line=False, caps=False),
                     ls_def,
                     ls_choices),
            Property("Width",
                     lambda x, ebar=ebar: outils.apply_to_ebar(
                         "set_linewidth", ebar, x, line=False),
                     bar.get_linewidth()[0],
                     None),
            Property("Color",
                     lambda x, ebar=ebar: outils.apply_to_ebar(
                         "set_color", ebar, x, line=False),
                     bar.get_color()[0],
                     None),
        ]


    return props, title


def _get_text_properties(txt):
    """ Properties for all text objects. """
    title = "Edit Text '{}'".format(txt.get_text())
    props = [
        Property("Label",  # needs to be called 'label' for fedit to know it's just text.
                 txt.set_text, txt.get_text(),
                 None,),
        "Style:",
        Property("Font",
                 lambda x, t=txt: _set_font(x, t), 
                 (txt.get_fontname(), int(txt.get_fontsize()),
                  "italic" == txt.get_fontstyle(), "normal" != txt.get_fontweight()
                  ),
                 None),
        Property("Color",
                 txt.set_color, txt.get_color(),
                 None),
        Property("Visible",
                 txt.set_visible, bool(txt.get_visible()),
                 None),
        "Position:",
        Property("X-Pos",
                 txt.set_x, txt.get_position()[0],
                 None),
        Property("Y-Pos",
                 txt.set_y, txt.get_position()[1],
                 None),
        Property("Z-Order",
                 txt.set_zorder, txt.get_zorder(),
                 None),
    ]

    return props, title


def _get_legend_properties(leg):
    title = "Edit Legend"
    props = [
        Property("Visible",
                 leg.set_visible, bool(leg.get_visible()),
                 None),
        Property("Draggable",
                 leg.draggable, leg._draggable is not None,
                 None),
        Property("Position",
                 lambda x, leg=leg: leg._set_loc(eval(x)),  # backconversion from string
                 str(leg._get_loc()),
                 None),
        Property("Z-Order",
                 leg.set_zorder, leg.get_zorder(),
                 None),
    ]
    return props, title


# Private Functions ############################################################


def _apply_data_to_properties(props, data, parent):
    """ Set the new data at the properties.

    data and props need to be in the same order, apart from the Separators in props.
    """
    props = [ppty for ppty in props if isinstance(ppty, Property)]  # remove separators
    for ppty, new_data in zip(props, data):
        ppty.set_value(new_data)

    # try to update the figure
    try:
        parent.update_figure()
    except AttributeError:
        pass


def _get_plane_from_axis(axis):
    if isinstance(axis, mpl.axis.XAxis):
        pos_choices = ["bottom", "top"]
        plane = 'x'
    else:
        pos_choices = ["left", "right"]
        plane = "y"
    return plane, pos_choices


def _set_font(font_tupel, txt):
    """ Helper to set font on text """
    txt.set_fontname(font_tupel[0])
    txt.set_fontsize(font_tupel[1])
    if font_tupel[2]:
        txt.set_fontstyle('italic')
    else:
        txt.set_fontstyle('normal')
    if font_tupel[3]:
        txt.set_fontweight('bold')
    else:
        txt.set_fontweight('normal')


def _set_axes_colors(color, axes):
    pass


def _set_axis_colors(color, axis):
    _, pos_choices = _get_plane_from_axis(axis)

    for spines in [axis.axes.spines[pos] for pos in pos_choices]:
        spines.set_color(color)

    ticks = axis.get_ticklines() + axis.get_ticklabels()
    for tick in ticks:
        tick.set_color(color)

    try:
        axis.label.set_color(color)
    except AttributeError:
        pass


def _show_ticks(state, axis, which):
    if which == "major":
        items = axis.get_major_ticks() + axis.get_majorticklabels()
    else:
        items = axis.get_minor_ticks() + axis.get_minorticklabels()

    for item in items:
        item.set_visible(state)



# Property Class ###############################################################


class Property(object):
    """ Object to automatically create the form-fields and apply their return value."""
    def __init__(self, title, setter, default, choices):
        self.setter = setter  # setter function
        self.title = title
        self.default = default
        self.choices = choices

    def get_form_line(self):
        default = self.default
        try:
            default = [default] + self.choices
        except TypeError:
            pass
        return self.title, default

    def set_value(self, value):
        self.setter(value)
