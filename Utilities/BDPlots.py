import traceback
import matplotlib
import matplotlib.pyplot as plt


class BDPlots:

    PREFIX = 'plots_handler__'

    def __init__(self, subplots_settings, plotter, logger):

        # Set logger
        self.logger = logger

        # Set plotter - the class that implements the plotting functions
        self.plotter = plotter

        # We initiate the subplots view to show them all (e.g. 0).
        # Once we zoom in, this will change to the relevant plot number
        self.subplots_view = 0

        # Subplots map
        self.subplots = {}
        self.sorted_subplots = []
        self.header_ax = None  # This is the axis we will "glue" the experiment header to

        # The subplots grid shape
        self.subplots_shape = subplots_settings['grid_shape']

        # Iterate over all subplots definitions and create a subplots map
        for subplot in subplots_settings['subplots']:

            # If the plot is not in the display list, skip creating an axis for it and mark it
            subplot['display'] = subplot["id"] in subplots_settings['display']
            if not subplot['display']:
                continue

            self.subplots[subplot["id"]] = subplot

        # Take the axis which is first in the list (e.g. Top-Left)
        self.subplots_header = self.get_subplot_by_index(0)

        # Initialize matplotlib
        self._prepare_matplotlib()

        pass

    def get_subplot_by_index(self, index):

        subplots_arr = []
        for subplot in self.subplots.values():
            if 'twinx' not in subplot.keys():
                subplots_arr.append(subplot)

        # Sort
        subplots_arr.sort(key=lambda sp: f'{sp["locy"]}{sp["locx"]}')

        return subplots_arr[index]

    def _prepare_matplotlib(self):

        # Ensure the backend we need
        matplotlib_version = matplotlib.get_backend()
        self.logger.info(f'Matplotlib backend: {matplotlib_version}')
        #matplotlib.use("Qt5Agg")

        pass

    def create_figures(self, new_figure=True):

        # Create the figure
        if new_figure:
            self.fig = plt.figure()

        # Create the subplots
        for subplot in self.subplots.values():
            if 'twinx' not in subplot.keys():
                subplot["ax"] = plt.subplot2grid(shape=self.subplots_shape,
                                           loc=(subplot['locy'], subplot['locx']),
                                           colspan=subplot['colspan'], rowspan=subplot['rowspan'])
            else:
                subplot["ax"] = self.subplots[subplot['twinx']]["ax"].twinx()

        pass

    def get_figure(self):
        return self.fig

    def should_display(self, subplot):
        """
        Returns True/False based on whether this subplot should be displayed
        """
        pass

    def set_figure_title(self, title):
        self.fig.canvas.manager.set_window_title(title)

    def plot_figures(self):

        if self.subplots_view > 0:
            subplots_values = [self.get_subplot_by_index(self.subplots_view-1)]
        else:
            subplots_values = self.subplots.values()

        # Iterate over all required figures and invoke their plotting function
        for subplot in subplots_values:

            # Check that we need to display this plot
            if not subplot["display"]:
                continue

            try:
                ax = subplot["ax"]
                ax.clear()

                def func_not_found():  # just in case we dont have the function
                    self.warn(f'Could not find plot function to invoke. Skipping')

                # Set subplot title
                if "title" in subplot.keys():
                    ax.set_title(subplot["title"], fontweight="bold")

                func_name = BDPlots.PREFIX + subplot['func']
                func = getattr(self.plotter, func_name, func_not_found)
                func(subplot)

                # Set subplot title and legend
                if "legend_loc" in subplot.keys():
                    ax.legend(loc=subplot["legend_loc"])

            except Exception as err:
                tb = traceback.format_exc()
                self.logger.warn(f'Problem with subplot {subplot["id"]}. Could not display it.\n{tb}')

        # Print the main experiment header (this is "glued" to the plot at [0,0]
        plot_header_func = getattr(self.plotter, BDPlots.PREFIX + "header", func_not_found)
        plot_header_func(self.subplots_header["ax"])

        # Plot the left sidebar
        plot_left_sidebar_func = getattr(self.plotter, BDPlots.PREFIX + "left_sidebar", func_not_found)
        plot_left_sidebar_func(self.subplots_header["ax"])

        # plt.tight_layout()
        plt.pause(0.2)

    def get_subplots_view(self):
        return self.subplots_view

    def set_subplots_view(self, new_subplot_view):
        """
        0 for all subplots
        [1, 2, ... n] - for a specific plot
        """

        # If there's no change, no need to do anything
        if new_subplot_view == self.subplots_view:
            return

        # If we're moving to a specific view, remove all subplots
        if self.subplots_view == 0 and new_subplot_view > 0:
            self.move_from_multiple_plots_to_single_plot(new_subplot_view)
        elif self.subplots_view > 0 and new_subplot_view == 0:
            self.move_from_single_plot_to_multiple_plots()
        elif self.subplots_view > 0 and new_subplot_view > 0:
            self.change_single_plot_to_other_single_plot(new_subplot_view)
        else:
            pass

        # Set the new mode to indicate we're showing this single subplot
        self.subplots_view = new_subplot_view

    def move_from_multiple_plots_to_single_plot(self, new_subplot_view):

        # Get the specific subplot we need to stay with
        subplot = self.get_subplot_by_index(new_subplot_view - 1)

        # Remove all other plots
        allaxes = self.fig.get_axes()
        for ax in allaxes:
            if ax.get_title() != subplot['title']:
                ax.remove()

        # Create new axis of the specific plot on the entire grid
        ax = plt.subplot2grid(shape=self.subplots_shape,
                              loc=(0, 0), colspan=self.subplots_shape[1], rowspan=self.subplots_shape[0])

        self.subplots[subplot['id']]['ax'] = ax

        # Set the subplot_header to this plot we zoomed into
        self.subplots_header = subplot

        pass

    def move_from_single_plot_to_multiple_plots(self):

        # Recreate all figures based on their original setting
        self.create_figures(new_figure=False)

        # Take the axis which is first in the list (e.g. Top-Left)
        self.subplots_header = self.get_subplot_by_index(0)

        pass

    def change_single_plot_to_other_single_plot(self, new_subplot_view):

        current_subplot = self.get_subplot_by_index(self.subplots_view - 1)
        new_subplot = self.get_subplot_by_index(new_subplot_view - 1)

        self.subplots[new_subplot['id']]['ax'] = current_subplot['ax']

        # Set the subplot_header to this plot we zoomed into
        self.subplots_header = new_subplot

        pass