import os
import json
import sys


class BDMenu:

    def __init__(self, caller, menu_file, menu_json=None):

        self.caller = caller

        if menu_file:
            self._load_menu(menu_file)
        elif menu_json:
            self.menu_data = menu_json

        pass

    def _load_menu(self, menu_file):

        try:
            # Opening JSON file
            f = open(menu_file)

            # returns JSON object as
            # a dictionary
            self.menu_data = json.load(f)

            # Closing file
            f.close()
        except Exception as err:
            print(err)

    def func_not_found(): # just in case we dont have the function
        print ('No Function Found!')

    def display(self):

        while True:
            self._display()
            if "infinite" not in self.menu_data or self.menu_data['infinite']==False:
                break
        print('Bye!')

    def _display(self):

        # TODO: check validity of json menu items

        call_map = {}

        # Iterating through the json menu items
        for menu_item in self.menu_data['menu_items']:
            if 'spacer' in menu_item:
                print()
            else:
                print(f'{menu_item["order"]}) {menu_item["display"]}')
                call_map[menu_item["order"]] = menu_item["func"]

        # Wait for input
        selection = input(">> ")
        selection = int(selection)

        # Get the function name and invoke it
        func_name = call_map[selection]

        if func_name.lower() == 'exit':
            sys.exit("User exited")

        func = getattr(self.caller, func_name, self.func_not_found)
        result = func()

        return result