import sys
from PyQt4 import QtGui
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt


class SbSGuiMatchResultView(QtGui.QWidget):

    def __init__(self, controller, variables_for_beam, variables_common, parent=None):
        super(SbSGuiMatchResultView, self).__init__(parent)

        self._latest_annotation = None
        self._disabled_constraints = {1: [], 2: []}
        self._elements_positions = controller.get_elements_positions()
        self._controller = controller
        self._variables_for_beam = variables_for_beam
        self._variables_common = variables_common
        self._build_gui()

    def _build_gui(self):
        main_layout = QtGui.QHBoxLayout()

        beam1_layout = QtGui.QVBoxLayout()
        beam1_frame = _BorderedGroupBox("Beam 1")
        beam1_frame.setLayout(beam1_layout)

        self._beam1_upper_figure = plt.figure()
        beam1_layout.addLayout(self._get_new_canvas_layout(self._beam1_upper_figure, 1))
        self._beam1_lower_figure = plt.figure()
        beam1_layout.addLayout(self._get_new_canvas_layout(self._beam1_lower_figure, 1))
        main_layout.addWidget(beam1_frame)

        beam2_layout = QtGui.QVBoxLayout()
        beam2_frame = _BorderedGroupBox("Beam 2")
        beam2_frame.setLayout(beam2_layout)

        self._beam2_upper_figure = plt.figure()
        beam2_layout.addLayout(self._get_new_canvas_layout(self._beam2_upper_figure, 2))
        self._beam2_lower_figure = plt.figure()
        beam2_layout.addLayout(self._get_new_canvas_layout(self._beam2_lower_figure, 2))
        main_layout.addWidget(beam2_frame)

        variables_layout = QtGui.QVBoxLayout()
        variables_frame = _BorderedGroupBox("Variables")
        variables_frame.setLayout(variables_layout)

        self._beam1_vars_layout = QtGui.QVBoxLayout()
        if not len(self._variables_for_beam[1]) == 0:
            beam1_widget = QtGui.QWidget()
            beam1_widget.setLayout(self._beam1_vars_layout)
            beam1_scroll = QtGui.QScrollArea()
            beam1_scroll.setWidgetResizable(True)
            beam1_scroll.setWidget(beam1_widget)
            beam1_vars_frame = _BorderedGroupBox("Beam 1")
            beam1_vars_frame_layout = QtGui.QHBoxLayout()
            beam1_vars_frame_layout.addWidget(beam1_scroll)
            beam1_vars_frame.setLayout(beam1_vars_frame_layout)
            variables_layout.addWidget(beam1_vars_frame)
            for variable in self._variables_for_beam[1]:
                self._beam1_vars_layout.addWidget(SbSGuiMatchResultView.CustomCheckBox(variable))

        self._beam2_vars_layout = QtGui.QVBoxLayout()
        if not len(self._variables_for_beam[2]) == 0:
            beam2_widget = QtGui.QWidget()
            beam2_widget.setLayout(self._beam2_vars_layout)
            beam2_scroll = QtGui.QScrollArea()
            beam2_scroll.setWidgetResizable(True)
            beam2_scroll.setWidget(beam2_widget)
            beam2_vars_frame = _BorderedGroupBox("Beam 2")
            beam2_vars_frame_layout = QtGui.QHBoxLayout()
            beam2_vars_frame_layout.addWidget(beam2_scroll)
            beam2_vars_frame.setLayout(beam2_vars_frame_layout)
            variables_layout.addWidget(beam2_vars_frame)
            for variable in self._variables_for_beam[2]:
                self._beam2_vars_layout.addWidget(SbSGuiMatchResultView.CustomCheckBox(variable))

        self._common_vars_layout = QtGui.QVBoxLayout()
        if not len(self._variables_common) == 0:
            common_widget = QtGui.QWidget()
            common_widget.setLayout(self._common_vars_layout)
            common_scroll = QtGui.QScrollArea()
            common_scroll.setWidgetResizable(True)
            common_scroll.setWidget(common_widget)
            common_vars_frame = _BorderedGroupBox("Common")
            common_vars_frame_layout = QtGui.QHBoxLayout()
            common_vars_frame_layout.addWidget(common_scroll)
            common_vars_frame.setLayout(common_vars_frame_layout)
            variables_layout.addWidget(common_vars_frame)
            for variable in self._variables_common:
                self._common_vars_layout.addWidget(SbSGuiMatchResultView.CustomCheckBox(variable))

        select_all_checkbox = QtGui.QCheckBox("Toggle select all")
        select_all_checkbox.stateChanged.connect(self._toogle_select_all)
        variables_layout.addWidget(select_all_checkbox)

        main_layout.addWidget(variables_frame)

        self.setLayout(main_layout)

    def _get_new_canvas_layout(self, figure, beam):
        layout = QtGui.QVBoxLayout()
        canvas = FigureCanvas(figure)
        canvas.mpl_connect(
            'motion_notify_event',
            lambda event: self._mouse_moved_on_figure(event, figure, beam)
        )
        canvas.mpl_connect(
            'button_press_event',
            lambda event: self._mouse_clicked_on_figure(event, figure, beam)
        )
        toolbar = SbSGuiMatchResultView.CustomNavigationBar(canvas, figure, self)
        layout.addWidget(toolbar)
        layout.addWidget(canvas)
        return layout

    def get_figures(self):
        return ((self._beam1_upper_figure, self._beam1_lower_figure),
                (self._beam2_upper_figure, self._beam2_lower_figure))

    def get_unselected_variables(self):
        unselected_vars = []

        def add_text_to_list(checkbox):
            if not checkbox.isChecked():
                unselected_vars.append(str(checkbox.text()))

        self._loop_through_checkboxes(add_text_to_list)
        return unselected_vars

    def get_disabled_constraints(self):
        disabled_constraints = []
        for constraints_list in self._disabled_constraints.values():
            disabled_constraints += constraints_list
        return disabled_constraints

    def _toogle_select_all(self, state):
        checked = bool(state)

        def toggle_checkbox(checkbox):
            checkbox.setChecked(checked)

        self._loop_through_checkboxes(toggle_checkbox)

    def update_variables(self, name_value_dict):
        def update(checkbox):
            checkbox.set_variable_value_from_dict(name_value_dict)

        self._loop_through_checkboxes(update)

    def _loop_through_checkboxes(self, function):
        for layout in [self._beam1_vars_layout,
                       self._beam2_vars_layout,
                       self._common_vars_layout]:
            for index in range(layout.count()):
                checkbox = layout.itemAt(index).widget()
                if issubclass(checkbox.__class__, QtGui.QCheckBox):
                    function(checkbox)

    def _mouse_moved_on_figure(self, event, figure, beam):
        for axes in figure.axes:
            del axes.texts[:]
        annotation_was_none = True
        if self._latest_annotation is not None:
            self._latest_annotation = None
            annotation_was_none = False
        axes, element_name, element_position = self.get_element_within_range(
            figure, event, beam
        )
        if element_name is not None and element_position is not None:
            x, y = element_position
            new_text = (element_name + "\n" +
                        "S = " + str(x))
            self._latest_annotation = axes.text(
                x, y,
                new_text,
                bbox=SbSGuiMatchResultController.BOX_STYLE
            )
        # TODO: Only redraw on element transition (no repeat same element)
        have_to_redraw = (
            (self._latest_annotation is not None) or
            (self._latest_annotation is None and not annotation_was_none)
        )
        if have_to_redraw:
            self._redraw_figure(beam, axes)

    def _mouse_clicked_on_figure(self, event, figure, beam):
        axes, selected_point = SbSGuiMatchResultView._get_point_within_range(
            figure,
            event,
        )
        if axes is not None and selected_point is not None:
            x, _ = selected_point
            element_name = self._elements_positions[beam][x]
            if element_name in self._disabled_constraints[beam]:
                self._disabled_constraints[beam].remove(element_name)
            else:
                self._disabled_constraints[beam].append(element_name)
        self._redraw_figure(beam, axes)

    def get_element_within_range(self, figure, event, beam):
        axes, selected_point = SbSGuiMatchResultView._get_point_within_range(
            figure,
            event,
        )
        if axes is not None and selected_point is not None:
            elements_positions = self._elements_positions[beam]
            x, y = selected_point
            element_name = elements_positions[x]
            element_position = (x, y)
            return axes, element_name, element_position
        return None, None, None

    def _get_points_for_element(self, element_name, beam, axes):
        points = []
        if len(axes.get_lines()) < 3:
            return []
        line = axes.get_lines()[2]
        xydata_in_plot = line.get_xydata()
        found_data_point = None
        for data_point in xydata_in_plot:
            x, _ = data_point
            if self._elements_positions[beam][x] == element_name:
                found_data_point = data_point
                break
        if found_data_point is not None:
            points.append(data_point)
        return points

    def _redraw_figure(self, beam, selected_axes):
        upper_figure = getattr(self, "_beam" + str(beam) + "_upper_figure")
        lower_figure = getattr(self, "_beam" + str(beam) + "_lower_figure")
        for figure in [upper_figure, lower_figure]:
            for axes in figure.axes:
                del axes.texts[:]
                for element_name in self._disabled_constraints[beam]:
                    points = self._get_points_for_element(element_name, beam, axes)
                    for point in points:
                        x, y = point
                        axes.text(x, y, "X",
                        horizontalalignment='center',
                        verticalalignment='center',
                        fontsize=15, color='red')
            if selected_axes is not None:
                selected_axes.texts.append(self._latest_annotation)
            figure.canvas.draw()

    @staticmethod
    def _get_point_within_range(figure, event):
        x_plot, y_plot = event.xdata, event.ydata
        if x_plot is not None and y_plot is not None:
            x_fig, y_fig = event.x, event.y
            for axes in figure.axes:
                for line in axes.get_lines():
                    xydata_in_plot = line.get_xydata()
                    min_distance2 = sys.float_info.max
                    selected_point_plot = None
                    for data_point in xydata_in_plot:
                        fig_point_x, fig_point_y = axes.transData.transform(data_point)
                        distance2 = (fig_point_x - x_fig) ** 2 + (fig_point_y - y_fig) ** 2
                        if distance2 < min_distance2:
                            min_distance2 = distance2
                            selected_point_plot = data_point
                    if min_distance2 < SbSGuiMatchResultController.DISTANCE_THRESHOLD2:
                        return axes, selected_point_plot
        return None, None

    class CustomNavigationBar(NavigationToolbar):
        LEGEND_LOCATIONS = [1, 2, 3, 4, None]

        def __init__(self, canvas, figure, parent=None):
            self.toolitems = list(self.toolitems)
            self.toolitems.append((
                "Move legend", "Move or hide the legend", "hand", "move_legend"
            ))
            super(SbSGuiMatchResultView.CustomNavigationBar, self).__init__(canvas, parent)

            self._figure = figure
            self._legend_location_index = 0

        def move_legend(self, *args):
            self._legend_location_index = (
                self._legend_location_index + 1
            ) % len(SbSGuiMatchResultView.CustomNavigationBar.LEGEND_LOCATIONS)
            for axes in self._figure.axes:
                legend = axes.get_legend()
                loc = SbSGuiMatchResultView.CustomNavigationBar.LEGEND_LOCATIONS[
                    self._legend_location_index
                ]
                if loc is None:
                    legend.set_visible(False)
                else:
                    legend.set_visible(True)
                    legend._set_loc(loc)
            self.draw()

    class CustomCheckBox(QtGui.QCheckBox):
        STRENGH_COLOR_LIMIT = 1e-3
        BACKGROUND_COLOR_TEMPLATE = "QCheckBox { background-color: %(COLOR)s;}"
        DISABLED_CSS_COLOR = "rgb(224, 224, 224)"

        def __init__(self, text):
            super(SbSGuiMatchResultView.CustomCheckBox, self).__init__(text)
            self.setMouseTracking(True)
            self._tooltip_text = ""

        def mouseMoveEvent(self, event):
            super(SbSGuiMatchResultView.CustomCheckBox, self).mouseMoveEvent(event)
            QtGui.QToolTip.showText(event.globalPos(),
                                    self._tooltip_text,
                                    widget=self)

        def set_variable_value_from_dict(self, name_value_dict):
            StaticRef = SbSGuiMatchResultView.CustomCheckBox
            try:
                variable_value = name_value_dict[str(self.text())]
                tooltip_text = str(variable_value)
                color_string = StaticRef._get_css_color_for_value(variable_value)
            except KeyError:
                tooltip_text = "Disabled"
                color_string = StaticRef.DISABLED_CSS_COLOR
            self._tooltip_text = str(tooltip_text)
            style = StaticRef.BACKGROUND_COLOR_TEMPLATE % {
                "COLOR": color_string,
            }
            self.setStyleSheet(style)
            self.update()

        @staticmethod
        def _get_css_color_for_value(value):
            absolute_strength = abs(value)
            limit = SbSGuiMatchResultView.CustomCheckBox.STRENGH_COLOR_LIMIT
            if absolute_strength > limit:
                red = 255
            else:
                red = int((absolute_strength / limit) * 255)
            green = 255 - red
            blue = 0
            return "rgb(" + str(red) + ", " + str(green) + ", " + str(blue) + ")"


class _BorderedGroupBox(QtGui.QGroupBox):

    def __init__(self, label, parent=None):
        super(_BorderedGroupBox, self).__init__(label, parent)


class SbSGuiMatchResultController(object):

    DISTANCE_THRESHOLD2 = 10 ** 2
    BOX_STYLE = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)

    def __init__(self, variables_for_beam, variables_common, get_positions_function):
        self._elements_positions = {}
        for beam in [1, 2]:
            try:
                self._elements_positions[beam] = get_positions_function(beam)
            except KeyError:
                continue
        self._view = SbSGuiMatchResultView(self, variables_for_beam, variables_common)

    def get_view(self):
        return self._view

    def get_unselected_variables(self):
        return self._view.get_unselected_variables()

    def get_disabled_constraints(self):
        return self._view.get_disabled_constraints()

    def get_figures(self):
        return self._view.get_figures()

    def get_elements_positions(self):
        return self._elements_positions

    def update_variables(self, name_value_dict):
        self._view.update_variables(name_value_dict)


if __name__ == "__main__":
    print >> sys.stderr, "This module is meant to be imported."
    sys.exit(-1)
