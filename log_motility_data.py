#!/usr/bin/env python

import os
import re
import sys
from argparse import ArgumentParser
from io import StringIO

import pandas as pd

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
rex_traj = re.compile(r"[\s\S]*?(?=Parameters)")
#\n\d[\s\S]*?(?=TotalSperm)
rex = re.compile(r"(?:\w+ +){7}\w+\n(?:(?:\d+\.\d+ {4}){7}\d+\.\d+\s)+")
rex_comp = re.compile(r"\d\nTotalSperm: \d+\.\d+\n((?:.+\s)+)#{13}.+#{13}")


def get_args():
    parser = ArgumentParser()
    parser.add_argument("filename")
    parser.add_argument("ID1")
    parser.add_argument("ID2")
    parser.add_argument("ID3")
    parser.add_argument("ID4", nargs='?', default='-')
    parser.add_argument("-o", "--output", nargs=1, default='Motility_Results.txt')
    args = parser.parse_args()
    return args


def main():
    # Use a breakpoint in the code line below to debug your script.
    args = get_args()
    filename_input = args.filename
    ID1 = args.ID1
    ID2 = args.ID2
    ID3 = args.ID3
    ID4 = args.ID4
    filename_output = args.output
    filename_output_t = 'traj.txt'

    # Read file
    if not os.path.isfile(filename_input):
        print("Input filename is not valid")
        sys.exit(0)
    with open(filename_input, "r") as f:
        text = f.read()
    m = rex.search(text)
    traj = rex_traj.search(text)
    m_comp = rex_comp.search(text)
    if m is None or m_comp is None or traj is None:
        print("Data not logged: Input file incomplete")
        sys.exit(0)
    test_data = StringIO(m.group())
    test_data_t = StringIO(traj.group())
    df = pd.read_csv(test_data, delim_whitespace=True)
    df_t = pd.read_csv(test_data_t, delim_whitespace=True, header=None)
    header_flag = not os.path.isfile(filename_output)
    df.insert(0, "ID4", ID4)
    df.insert(0, "ID3", ID3)
    df.insert(0, "ID2", ID2)
    df.insert(0, "ID1", ID1)
    df_t.insert(0, "ID4", ID4)
    df_t.insert(0, "ID3", ID3)
    df_t.insert(0, "ID2", ID2)
    df_t.insert(0, "ID1", ID1)
    df.to_csv(filename_output, mode='a', header=header_flag, index=False)
    df_t.to_csv(filename_output_t, mode='a', header=None, index=False)
    string_comp = m_comp.group(1)
    list_comp = string_comp.split('\n')

    for i in range(int(len(list_comp) / 2)):
        header = list_comp[2 * i]
        values = list_comp[2 * i + 1]

        base_name, _ = os.path.splitext(filename_output)
        first_header = header.split()[0]
        if first_header == "TM":
            filename_comp = f"{base_name}_Avg.txt"
        elif first_header == "SLOW":
            filename_comp = f"{base_name}_Speed.txt"
        elif first_header == "MedianVCL":
            filename_comp = f"{base_name}_Median.txt"
        elif first_header == "SigmaVCL":
            filename_comp = f"{base_name}_Sigma.txt"
        else:
            breakpoint()
            raise ValueError(f"Unrecognized header: {first_header}")

        header_flag = not os.path.isfile(filename_comp)
        with open(filename_comp, "a") as f:
            if header_flag:
                f.write("ID1 ID2 ID3 ID4 " + header + "\n")
            f.write(ID1 + " " + ID2 + " " + ID3 + " " + ID4 + " " + values + "\n")


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
