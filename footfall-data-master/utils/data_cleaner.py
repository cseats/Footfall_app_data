import os
import pandas as pd
import re
import decimal
import json


def get_files(path):
    files = []
    for f in os.listdir(path):
        if f.startswith("data"):
            files.append(f)

    return files


def clean(path):
    files = get_files(path)
    for file in files:
        # data = pd.read_excel(path + '\\' + file, 'Sheet1')
        data = open(path + '\\' + file, 'r')

    data = data.read().split('\n')
    data.pop(0)

    print(len(data))

    accels = []

    for rec in data:
        rec = rec.split('\t')
        points = []
        if rec[0] != '' and rec[1] != '':
            if ';;' in rec[1]:
                rec[1] = rec[1].split(';;')
                for point in rec[1]:
                    point = point.split(',')
                    if len(point) > 2:
                        temp1 = point[0] + '.' + point[1]
                        temp2 = point[2] + '.' + point[3]
                        point = [temp1, temp2]
                    points.append({'time': float(point[0]), 'acceleration': float(point[1])})
                accels.append(points)
    print(len(accels))
    with open('clean_data.json', 'w') as outfile:
        json.dump(accels, outfile)

    return accels


def clean_from_database(data):
    print(len(data))

    accels = []

    for rec in data:
        points = []
        if rec['RecordedData']:
            if ';;' in rec['RecordedData']:
                rec['RecordedData'] = rec['RecordedData'].split(';;')
                for point in rec['RecordedData']:
                    point = point.split(',')
                    if len(point) == 4:
                        temp1 = point[0] + '.' + point[1]
                        temp2 = point[2] + '.' + point[3]
                        point = [temp1, temp2]
                    if len(point) == 2:
                        points.append({
                            'time': float(point[0]),
                            'acceleration_x': None,
                            'acceleration_y': None,
                            'acceleration_z': float(point[1])
                        })
                    elif len(point) > 2:
                        points.append({
                            'time': float(point[0]),
                            'acceleration_x': float(point[1]),
                            'acceleration_y': float(point[2]),
                            'acceleration_z': float(point[3])
                        })
                accels.append(points)
    print(len(accels))
    with open('clean_db_data.json', 'w') as outfile:
        json.dump(accels, outfile)

    return accels
