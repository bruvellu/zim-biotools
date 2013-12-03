# -*- coding: utf-8 -*-

# Bruno C. Vellutini <organelas@gmail.com>

'''Bioinformatics plugin for Zim. [WARNING, NOT FUNCTIONAL]'''

import gtk

from zim.plugins import PluginClass

ui_xml = '''
<ui>
    <menubar name='menubar'>
        <menu action='tools_menu'>
            <placeholder name='plugin_items'>
                <menuitem action='fetch_ncbi'/>
            </placeholder>
        </menu>
    </menubar>
</ui>
'''

ui_actions = (
    # name, stock id, label, accelerator, tooltip, readonly
    ('fetch_ncbi', 'gtk-ncbi', _('_Fetch sequence'), '<ctrl><alt>F', 'Fetch NCBI', True),  # T: menu item

)


class BioTools(PluginClass):

    plugin_info = {
        'name': _('BioTools'),  # T: plugin name
        'description': _('This plugin provides useful bioinformatics tools.'),  # T: plugin description
        'author': 'Bruno C. Vellutini',
    }

    def __init__(self, ui):
        PluginClass.__init__(self, ui)

        # Add menu items.
        if self.ui.ui_type == 'gtk':
            self.ui.add_actions(ui_actions, self)
            self.ui.add_ui(ui_xml, self)

    def fetch_ncbi(self, ref=None):
        '''Fetch aminoacid sequences using NCBI identifier.'''
        # Get selection, otherwise abort!
        buffer = self.ui.mainwindow.pageview.view.get_buffer()
        print buffer
