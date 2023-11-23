import os
import json
import sys
import numpy as np


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
                call_map[menu_item["order"]] = menu_item

        # Wait for input
        selection = input(">> ")
        selection = int(selection)

        # Sanity check
        if "action" not in call_map[selection]:
            self._error('Missing "action" section in menu definitions. Check your json file.')

        # Get the function name and invoke it
        func_name = call_map[selection]["action"]

        # If it is an exit function, then exit
        if func_name.lower() == 'exit':
            sys.exit("User exited")
            return

        args = {}
        if "args" in call_map[selection]:
            for arg in call_map[selection]["args"]:
                value_str = input(arg['name'])
                value = self._convert_input(value_str, arg['type'])
                args[arg['name']] = value

        # Invoke the method/action
        func = getattr(self.caller, func_name, self.func_not_found)

        if len(args) == 0:
            result = func()
        else:
            result = func(**args)

        return result

    def _convert_input(self, value_str, value_type):
        if 'array_of' in value_type:
            # Remove redundant spaces and split by commas
            value_str = value_str.replace(" ", "").split(',')
            # Get the type requested
            value_type = value_type.replace('array_of_', '')
            # Create array
            values = []
            for val in value_str:
                values.append(self._convert_input(val, value_type))
            return values
        elif value_type == 'int' or value_type == 'integer':
            return int(float(value_str))  # Python fails to convert float strings with this: int("3.1415"), so we need to convert to floar first and then int
        elif value_type == 'float':
            return float(value_str)
        elif value_type == 'string' or value_type == 'str':
            return str(value_str)
        else:
            self._error(f'Unknown value_type ({value_type}). Check your json.')
        pass
    def _error(self, message):
        print(message)
        sys.exit("User exited")
        pass