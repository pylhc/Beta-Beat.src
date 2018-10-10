# -*- coding: utf-8 -*-
#
# Copyright © 2009 Pierre Raybaut
# Licensed under the terms of the MIT License
# see the mpl licenses directory for a copy of the license
#
# jdilly: I do not like it.
#         It's very hard to modularize the way it is, but I did my best.
#         Should be rewritten from scratch.

"""Module that provides a GUI-based editor for matplotlib's figure options"""

from __future__ import (division, print_function,
                        unicode_literals)

import re

import matplotlib
from matplotlib import cm, colors as mcolors, image as mimage
import matplotlib.backends.qt_editor.formlayout as formlayout

import options_utils as outils
from gui_utils import get_icon

LINESTYLES = outils.get_linestyles()
DRAWSTYLES = outils.get_drawstyles()
MARKERS = outils.get_markers()
NOLEGEND = outils.NOLEGEND


def figure_edit(axes, parent=None):
    """Edit matplotlib figure options"""
    sep = (None, None)  # separator

    # Get / General
    # Cast to builtin floats as they have nicer reprs.
    xmin, xmax = map(float, axes.get_xlim())
    ymin, ymax = map(float, axes.get_ylim())
    general = [('Title', axes.get_title()),
               sep,
               (None, "<b>X-Axis</b>"),
               ('Left', xmin), ('Right', xmax),
               ('Label', axes.get_xlabel()),
               ('Scale', [axes.get_xscale(), 'linear', 'log', 'logit']),
               sep,
               (None, "<b>Y-Axis</b>"),
               ('Bottom', ymin), ('Top', ymax),
               ('Label', axes.get_ylabel()),
               ('Scale', [axes.get_yscale(), 'linear', 'log', 'logit']),
               sep,
               ('(Re-)Generate automatic legend', False),
               ]

    # Save the unit data
    xconverter = axes.xaxis.converter
    yconverter = axes.yaxis.converter
    xunits = axes.xaxis.get_units()
    yunits = axes.yaxis.get_units()

    # Sorting for default labels (_lineXXX, _imageXXX).
    def cmp_key(label):
        match = re.match(r"(_line|_image)(\d+)", label)
        if match:
            return match.group(1), int(match.group(2))
        else:
            return label, 0

    # Get / Curves
    linedict = {}
    for line in axes.get_lines():
        label = line.get_label()
        if label == NOLEGEND:
            continue
        linedict[label] = line

    for idx, ebar in enumerate([cont for cont in axes.containers if outils.is_errorbar(cont)]):
            label = ebar.get_label()
            if label == NOLEGEND:
                label = "unnamed{:d}".format(idx)
            linedict[label] = ebar

    curves = []

    curvelabels = sorted(linedict, key=cmp_key)
    for label in curvelabels:
        line = linedict[label]
        if outils.is_errorbar(line):
            line = line[0]

        color = mcolors.to_hex(
            mcolors.to_rgba(line.get_color(), line.get_alpha()),
            keep_alpha=True)
        ec = mcolors.to_hex(
            mcolors.to_rgba(line.get_markeredgecolor(), line.get_alpha()),
            keep_alpha=True)
        fc = mcolors.to_hex(
            mcolors.to_rgba(line.get_markerfacecolor(), line.get_alpha()),
            keep_alpha=True)
        ls_data = outils.prepare_formdata(LINESTYLES, line.get_linestyle())
        ds_data = outils.prepare_formdata(DRAWSTYLES, line.get_drawstyle())
        ms_data = outils.prepare_formdata(MARKERS, line.get_marker())
        curvedata = [
            ('Label', label),
            ('Visible', line.get_visible()),
            sep,
            (None, '<b>Line</b>'),
            ('Line style', [ls_data[0]] + ls_data[1]),
            ('Draw style', [ds_data[0]] + ds_data[1]),
            ('Width', line.get_linewidth()),
            ('Color (RGBA)', color),
            sep,
            (None, '<b>Marker</b>'),
            ('Style', [ms_data[0]] + ms_data[1]),
            ('Size', line.get_markersize()),
            ('Face color (RGBA)', fc),
            ('Edge color (RGBA)', ec),
        ]

        curves.append([curvedata, label, ""])
    # Is there a curve displayed?
    has_curve = bool(curves)

    # Get / Images
    imagedict = {}
    for image in axes.get_images():
        label = image.get_label()
        if label == NOLEGEND:
            continue
        imagedict[label] = image
    imagelabels = sorted(imagedict, key=cmp_key)
    images = []
    cmaps = [(cmap, name) for name, cmap in sorted(cm.cmap_d.items())]
    for label in imagelabels:
        image = imagedict[label]
        cmap = image.get_cmap()
        if cmap not in cm.cmap_d.values():
            cmaps = [(cmap, cmap.name)] + cmaps
        low, high = image.get_clim()
        imagedata = [
            ('Label', label),
            ('Colormap', [cmap.name] + cmaps),
            ('Min. value', low),
            ('Max. value', high),
            ('Interpolation',
             [image.get_interpolation()]
             + [(name, name) for name in sorted(mimage.interpolations_names)])]
        images.append([imagedata, label, ""])
    # Is there an image displayed?
    has_image = bool(images)

    datalist = [(general, "Axes", "")]
    if curves:
        datalist.append((curves, "Curves", ""))
    if images:
        datalist.append((images, "Images", ""))

    def apply_callback(data):
        """This function will be called to apply changes"""
        orig_xlim = axes.get_xlim()
        orig_ylim = axes.get_ylim()

        general = data.pop(0)
        curves = data.pop(0) if has_curve else []
        images = data.pop(0) if has_image else []
        if data:
            raise ValueError("Unexpected field")

        # Set / General
        (title, xmin, xmax, xlabel, xscale, ymin, ymax, ylabel, yscale,
         generate_legend) = general

        if axes.get_xscale() != xscale:
            axes.set_xscale(xscale)
        if axes.get_yscale() != yscale:
            axes.set_yscale(yscale)

        axes.set_title(title)
        axes.set_xlim(xmin, xmax)
        axes.set_xlabel(xlabel)
        axes.set_ylim(ymin, ymax)
        axes.set_ylabel(ylabel)

        # Restore the unit data
        axes.xaxis.converter = xconverter
        axes.yaxis.converter = yconverter
        axes.xaxis.set_units(xunits)
        axes.yaxis.set_units(yunits)
        axes.xaxis._update_axisinfo()
        axes.yaxis._update_axisinfo()

        # Set / Curves
        for index, curve in enumerate(curves):
            (label, visible, linestyle, drawstyle, linewidth, color, marker, markersize,
             markerfacecolor, markeredgecolor) = curve

            line = linedict[curvelabels[index]]
            hbar = None
            vbar = None
            ebar = None
            line.set_label(label)

            if outils.is_errorbar(line):
                ebar = line
                hbar = ebar[1]
                vbar = ebar[2]
                line = ebar[0]

            line.set_linestyle(linestyle)
            line.set_drawstyle(drawstyle)
            line.set_linewidth(linewidth)
            rgba = mcolors.to_rgba(color)
            line.set_alpha(None)
            line.set_color(rgba)
            line.set_visible(visible)
            for bar in [hbar, vbar]:
                if bar is not None:
                    for b in bar:
                        b.set_linewidth(linewidth)
                        b.set_color(rgba)
                        b.set_visible(visible)

            if marker is not 'none':
                line.set_marker(marker)
                line.set_markersize(markersize)
                line.set_markerfacecolor(markerfacecolor)
                line.set_markeredgecolor(markeredgecolor)

            if not visible:
                if ebar:
                    ebar.set_label(NOLEGEND)
                else:
                    line.set_label(NOLEGEND)

        # Set / Images
        for index, image_settings in enumerate(images):
            image = imagedict[imagelabels[index]]
            label, cmap, low, high, interpolation = image_settings
            image.set_label(label)
            image.set_cmap(cm.get_cmap(cmap))
            image.set_clim(*sorted([low, high]))
            image.set_interpolation(interpolation)

        # re-generate legend, if checkbox is checked
        if generate_legend:
            outils.regenerate_legend(axes)

        # Redraw
        figure = axes.get_figure()
        figure.canvas.draw()
        if not (axes.get_xlim() == orig_xlim and axes.get_ylim() == orig_ylim):
            figure.canvas.toolbar.push_current()

    data = formlayout.fedit(datalist, title="Figure options", parent=parent,
                            icon=get_icon('qt4_editor_options.svg'),
                            apply=apply_callback)
    if data is not None:
        apply_callback(data)

