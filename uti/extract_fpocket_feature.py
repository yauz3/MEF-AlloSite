# 23/02/2024
# Author: Sadettin Y. Ugurlu

import os


def extract_pocket_number(pocket):
    """
    The def search pocket id.
    :param pocket: Pocket NUMBER
    :return: POCKET NUMBER or None if a number is not found.
    """
    # find the pocket file path
    file_name = os.path.basename(pocket)
    # find the "pocket" in the file
    pocket_index = file_name.find('pocket')
    if pocket_index != -1:
        # take the part next to the "pocket"
        pocket_str = file_name[pocket_index + len('pocket'):]
        # find the number here
        pocket_num = ''
        for char in pocket_str:
            if char.isdigit():
                pocket_num += char
            else:
                break
        # if you find a number, return it
        if pocket_num:
            return int(pocket_num)
    # if there is no number, return None
    return None


def parse_pocket_values(file_path):
    """
    Function that converts data for each pocket from a text file to a dictionary with values.

    Args:
    file_path: Path of the text file.

    Returns:
    Pocket data contained in a dictionary or None if the file is not found.
    """
    pocket_data = {}
    with open(file_path, 'r') as file:
        lines = file.readlines()

        current_pocket = None
        for line in lines:
            stripped_line = line.strip()
            if not stripped_line:  # Empty line
                continue

            if stripped_line.startswith("Pocket"):
                current_pocket = stripped_line.split()[1]  # current pocket
                pocket_data[current_pocket] = {}
            else:
                key, value = stripped_line.split(':')
                pocket_data[current_pocket][key.strip()] = value.strip()

    return pocket_data

