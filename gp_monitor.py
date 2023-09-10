import gi

gi.require_version("Gtk", "3.0")

from matplotlib.pyplot import inferno
from time import time
from pypgtable import table
from logging import DEBUG, INFO, WARN, ERROR, FATAL, NullHandler, getLogger
from graph_tool import Graph, Vertex
from graph_tool.util import find_vertex
from graph_tool.draw import (
    graph_draw,
    interactive_window,
    sfdp_layout,
    GraphWindow,
    GraphWidget,
)
from graph_tool.topology import label_components
from matplotlib.cm import viridis, plasma, inferno, magma, cividis
from matplotlib.cm import (
    Greys as greys,
    Purples as purples,
    Blues as blues,
    Greens as greens,
    Oranges as oranges,
    Reds as reds,
)
from matplotlib.cm import cool, hot, gist_heat
from matplotlib.cm import (
    gnuplot,
    gnuplot2,
    CMRmap as cmrmap,
    rainbow,
    jet,
    nipy_spectral,
)
from gi.repository import Gtk, Gdk, GdkPixbuf, GObject, GLib, Gio

_logger = getLogger(__name__)
_logger.addHandler(NullHandler())


# The gene pool default config
# Tree structure
_LEL = "ancestor_a_ref"
_REL = "ancestor_b_ref"
_NL = "ref"
_PTR_MAP = {_LEL: _NL, _REL: _NL}
_DEFAULT_GP_CONFIG = {
    "database": {
        "dbname": "test_db",
    },
    "ptr_map": _PTR_MAP,
    "table": "gene_pool",
}
_COLUMNS = (_NL, _LEL, _REL, "survivability")

_OVER_MAX = 1 << 64
_MASK = _OVER_MAX - 1
ref_str = lambda x: f"{(_OVER_MAX + x) & _MASK:016x}"

print("Loading Gene Pool...")
gp = table(_DEFAULT_GP_CONFIG)
gcs = tuple(gp.recursive_select("WHERE {population_uid}={pid}", literals={"pid": 1}))
# for gc in gcs: print(gc)

print(f"Creating population of {len(gcs)} GC's in graph...")
g = Graph()
vertices = g.add_vertex(len(gcs))
pin = g.new_vertex_property("bool", val=False)
refstr = g.new_vertex_property("string", vals=(ref_str(gc[_NL]) for gc in gcs))
ref = g.new_vertex_property("int64_t", vals=(gc[_NL] for gc in gcs))
s = g.new_vertex_property("float", vals=(gc["survivability"] for gc in gcs))
for gc in gcs:
    v = find_vertex(g, ref, gc[_NL])[0]
    if gc[_LEL] is not None:
        c = find_vertex(g, ref, gc[_LEL])
        if c:
            g.add_edge(c[0], v)
        else:
            nv = g.add_vertex()
            refstr[nv] = ref_str(gc[_LEL])
            ref[nv] = gc[_LEL]
            s[nv] = 0.0
            g.add_edge(nv, v)

    if gc[_REL] is not None:
        c = find_vertex(g, ref, gc[_REL])
        if c:
            g.add_edge(c[0], v)
        else:
            nv = g.add_vertex()
            refstr[nv] = ref_str(gc[_REL])
            ref[nv] = gc[_REL]
            s[nv] = 0.0
            g.add_edge(nv, v)

print("Calculating initial layout...", end="", flush=True)
start = time()
count = 0
pos = sfdp_layout(g)
print(f"{time()-start:.2f}s")
print("Creating graph widget...")


class egp_GraphWidget(GraphWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.single_gc = False
        self.set_hexpand(True)
        self.set_vexpand(True)
        self.set_halign(Gtk.Align.FILL)
        self.set_valign(Gtk.Align.FILL)

    """Add a single left click event"""

    def button_release_event(self, widget, event):
        if event.button == 1 and not self.moved_picked:
            if not self.single_gc and isinstance(self.picked, Vertex):
                components, _ = label_components(self.g, directed=False)
                component = components[self.picked]
                view = g.new_vertex_property("bool", val=False)
                for v in find_vertex(self.g, components, component):
                    view[v] = True
                self.single_gc = True
                self.g.set_vertex_filter(view)
                self.regenerate_surface(reset=True)
                self.queue_draw()
        super().button_release_event(widget, event)

    def motion_notify_event(self, widget, event):
        state = event.window.get_pointer()[3] if event.is_hint else event.state
        if not (
            state & Gdk.ModifierType.BUTTON1_MASK
            and not state & Gdk.ModifierType.CONTROL_MASK
            and not state & Gdk.ModifierType.SHIFT_MASK
        ):
            super().motion_notify_event(widget, event)


# graph_draw(g, vcmap=inferno, vertex_fill_color=s, bg_color=(1,1,1,1), output_size=(1000,1000), output='gcs.png')
graph = egp_GraphWidget(
    g, fit_view=True, pos=pos, vcmap=jet, vertex_fill_color=s, bg_color=(1, 1, 1, 1)
)

import sys

# This would typically be its own file
MENU_XML = """
<?xml version="1.0" encoding="UTF-8"?>
<interface>
  <menu id="app-menu">
    <section>
      <item>
        <attribute name="action">app.about</attribute>
        <attribute name="label" translatable="yes">_About</attribute>
      </item>
      <item>
        <attribute name="action">app.quit</attribute>
        <attribute name="label" translatable="yes">_Quit</attribute>
        <attribute name="accel">&lt;Primary&gt;q</attribute>
      </item>
    </section>
  </menu>
</interface>
"""

CM_XML = """
<?xml version="1.0" encoding="UTF-8"?>
<interface>
  <menu id="cm-menu">
    <submenu>
        <attribute name="label" translatable="yes">Colour Mappings</attribute>
        <section>
        <item>
            <attribute name="action">win.colour_mapping</attribute>
            <attribute name="target">Viridis</attribute>
            <attribute name="label" translatable="yes">Viridis</attribute>
        </item>
        <item>
            <attribute name="action">win.colour_mapping</attribute>
            <attribute name="target">Inferno</attribute>
            <attribute name="label" translatable="yes">Inferno</attribute>
        </item>
        <item>
            <attribute name="action">win.colour_mapping</attribute>
            <attribute name="target">Plasma</attribute>
            <attribute name="label" translatable="yes">Plasma</attribute>
        </item>
        <item>
            <attribute name="action">win.colour_mapping</attribute>
            <attribute name="target">Magma</attribute>
            <attribute name="label" translatable="yes">Magma</attribute>
        </item>
        <item>
            <attribute name="action">win.colour_mapping</attribute>
            <attribute name="target">Cividis</attribute>
            <attribute name="label" translatable="yes">Cividis</attribute>
        </item>
        <item>
            <attribute name="action">win.colour_mapping</attribute>
            <attribute name="target">Cool</attribute>
            <attribute name="label" translatable="yes">Cool</attribute>
        </item>
        <item>
            <attribute name="action">win.colour_mapping</attribute>
            <attribute name="target">Hot</attribute>
            <attribute name="label" translatable="yes">Hot</attribute>
        </item>
        <item>
            <attribute name="action">win.colour_mapping</attribute>
            <attribute name="target">Gist_heat</attribute>
            <attribute name="label" translatable="yes">Gist_heat</attribute>
        </item>
        <item>
            <attribute name="action">win.colour_mapping</attribute>
            <attribute name="target">Gnuplot</attribute>
            <attribute name="label" translatable="yes">Gnuplot</attribute>
        </item>
        <item>
            <attribute name="action">win.colour_mapping</attribute>
            <attribute name="target">Gnuplot2</attribute>
            <attribute name="label" translatable="yes">Gnuplot2</attribute>
        </item>
        <item>
            <attribute name="action">win.colour_mapping</attribute>
            <attribute name="target">CMRmap</attribute>
            <attribute name="label" translatable="yes">CMRmap</attribute>
        </item>
        <item>
            <attribute name="action">win.colour_mapping</attribute>
            <attribute name="target">Rainbow</attribute>
            <attribute name="label" translatable="yes">Rainbow</attribute>
        </item>
        <item>
            <attribute name="action">win.colour_mapping</attribute>
            <attribute name="target">Jet</attribute>
            <attribute name="label" translatable="yes">Jet</attribute>
        </item>
        <item>
            <attribute name="action">win.colour_mapping</attribute>
            <attribute name="target">Nipy_spectral</attribute>
            <attribute name="label" translatable="yes">Nipy_spectral</attribute>
        </item>
        <item>
            <attribute name="action">win.colour_mapping</attribute>
            <attribute name="target">Greys</attribute>
            <attribute name="label" translatable="yes">Greys</attribute>
        </item>
        <item>
            <attribute name="action">win.colour_mapping</attribute>
            <attribute name="target">Purples</attribute>
            <attribute name="label" translatable="yes">Purples</attribute>
        </item>
        <item>
            <attribute name="action">win.colour_mapping</attribute>
            <attribute name="target">Blues</attribute>
            <attribute name="label" translatable="yes">Blues</attribute>
        </item>
        <item>
            <attribute name="action">win.colour_mapping</attribute>
            <attribute name="target">Greens</attribute>
            <attribute name="label" translatable="yes">Greens</attribute>
        </item>
        <item>
            <attribute name="action">win.colour_mapping</attribute>
            <attribute name="target">Oranges</attribute>
            <attribute name="label" translatable="yes">Oranges</attribute>
        </item>
        <item>
            <attribute name="action">win.colour_mapping</attribute>
            <attribute name="target">Reds</attribute>
            <attribute name="label" translatable="yes">Reds</attribute>
        </item>
        </section>
    </submenu>
  </menu>
</interface>
"""


class AppWindow(Gtk.ApplicationWindow):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.set_default_size(1200, 1000)

        colour_mapping_action = Gio.SimpleAction.new_stateful(
            "colour_mapping",
            GLib.VariantType.new("s"),
            GLib.Variant.new_string("Inferno"),
        )
        colour_mapping_action.connect("activate", self.set_colour_mapping)
        self.add_action(colour_mapping_action)

        self.grid = Gtk.Grid()

        self.label1 = Gtk.Label(margin=10)
        self.label1.set_markup(
            "(<b>R</b>)esize to window\n"
            "(<b>S</b>)pring-block layout regeneration\n"
            "(<b>Z</b>)oom to selection\n"
            "(<b>C</b>)lear family tree filter.\n\n"
            "Left click to view nearest family tree.\n"
            "Middle mouse drag to pan\n"
            "Scroll to zoom\n"
            "SHIFT + scroll to resize vertices &amp; edges\n"
            "SHIFT + left mouse drag to select vertices\n"
            "CRTL + scroll to rotate about pointer"
        )
        self.label1.set_justify(Gtk.Justification.LEFT)

        self.separator = Gtk.HSeparator()

        self.label2_text = "Reference\t\t: {0}\n" "Survivability\t: {1}\n"

        self.label2 = Gtk.Label(margin=10)
        self.label2.set_markup(self.label2_text.format("-", "-"))
        self.label2.set_justify(Gtk.Justification.LEFT)
        self.label2.set_halign(Gtk.Align.START)

        self.grid.attach(graph, 0, 0, 5, 9)
        self.grid.attach_next_to(self.label1, graph, Gtk.PositionType.RIGHT, 1, 4)
        self.grid.attach_next_to(
            self.separator, self.label1, Gtk.PositionType.BOTTOM, 1, 1
        )
        self.grid.attach_next_to(
            self.label2, self.separator, Gtk.PositionType.BOTTOM, 1, 4
        )

        self.add(self.grid)
        self.show_all()

        graph.key_press_user_callback = self.key_press_callback
        self.connect_after("motion-notify-event", self.motion_notify_event)
        self.previous_picked = None

    def set_colour_mapping(self, action, value):
        action.set_state(value)
        graph.kwargs["vcmap"] = eval(value.get_string().lower())
        graph.regenerate_surface(reset=True)
        graph.queue_draw()

    def key_press_callback(self, w, g, keyval, picked, pos, vprops, eprops):
        if keyval == ord("c"):
            graph.single_gc = False
            graph.g.set_vertex_filter(None)
            graph.regenerate_surface(reset=True)
            graph.queue_draw()

    """Add update panel information picked vertex."""

    def motion_notify_event(self, widget, event):
        if isinstance(graph.picked, Vertex):
            label2_text = self.label2_text.format(
                refstr[graph.picked], "{:0.3f}".format(s[graph.picked])
            )
            self.previous_picked = graph.picked
            self.label2.set_markup(label2_text)
        elif self.previous_picked is not None:
            label2_text = self.label2_text.format("-", "-")
            self.label2.set_markup(label2_text)
        # else did not change


class Application(Gtk.Application):
    def __init__(self, *args, **kwargs):
        super().__init__(
            *args,
            application_id="org.example.myapp",
            flags=Gio.ApplicationFlags.HANDLES_COMMAND_LINE,
            **kwargs,
        )
        self.window = None

        self.add_main_option(
            "test",
            ord("t"),
            GLib.OptionFlags.NONE,
            GLib.OptionArg.NONE,
            "Command line test",
            None,
        )

    def do_startup(self):
        Gtk.Application.do_startup(self)

        action = Gio.SimpleAction.new("about", None)
        action.connect("activate", self.on_about)
        self.add_action(action)

        action = Gio.SimpleAction.new("quit", None)
        action.connect("activate", self.on_quit)
        self.add_action(action)

        builder = Gtk.Builder.new_from_string(MENU_XML, -1)
        self.set_app_menu(builder.get_object("app-menu"))

        builder = Gtk.Builder.new_from_string(CM_XML, -1)
        self.set_menubar(builder.get_object("cm-menu"))

    def do_activate(self):
        # We only allow a single window and raise any existing ones
        if not self.window:
            # Windows are associated with the application
            # when the last one is closed the application shuts down
            self.window = AppWindow(application=self, title="Main Window")

        self.window.present()

    def do_command_line(self, command_line):
        options = command_line.get_options_dict()
        # convert GVariantDict -> GVariant -> dict
        options = options.end().unpack()

        if "test" in options:
            # This is printed on the main instance
            print("Test argument recieved: %s" % options["test"])

        self.activate()
        return 0

    def on_about(self, action, param):
        about_dialog = Gtk.AboutDialog(transient_for=self.window, modal=True)
        about_dialog.present()

    def on_quit(self, action, param):
        self.quit()


if __name__ == "__main__":
    print("Running application...")
    app = Application()
    app.run(sys.argv)
