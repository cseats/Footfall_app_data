import os
import json
import time
import pandas as pd
from vibration import processing
from utils import data_cleaner, database

dir_path = os.path.dirname(os.path.realpath(__file__))

# need to ensure that the ip addresses that are being used are approved by the firewall

def update_data():
    db = database.connect()
    query = """
    SELECT [Id]
          ,[RecordedData]
          ,[MeasuredPeakResponseFactor]
          ,[MeasuredNaturalFreq]
          ,[Metadata]
          ,[Location]
    FROM [dbo].[VibrationItems]
    """
    rows = database.run_pandas_query(db, query)
    rows.to_json('data\\db_out.json', orient='table')

    accels = data_cleaner.clean_from_database(rows)
    return accels


def get_good_data():
    db = database.connect()
    query = """
        SELECT top 10 [Id]
              ,[RecordedData]
              ,[MeasuredPeakResponseFactor]
              ,[MeasuredNaturalFreq]
              ,[Metadata]
              ,[Location]
        FROM [dbo].[VibrationItems]
        WHERE  [Location] IS NOT NULL AND [Metadata] IS NOT NULL
    """
    rows = database.run_pandas_query(db, query)
    print(rows)
    # rows.to_json('db_out.json', orient='table')

    # accels = data_cleaner.clean_from_database(rows)
    return


def count_db_entries():
    db = database.connect()
    query = """
        SELECT count(*) FROM [dbo].[VibrationItems]
        """
    count = database.run_query(db, query)
    print("Number of entries: %g" % count)

    return count


def fetch_and_clean():
    res = update_data()
    clean_data = data_cleaner.clean_from_database(res)

    return clean_data


def fetch_full_pipeline():
    res = fetch_and_clean()
    processing.processing(res)


def process_existing():
    # dirty_data = pd.read_json(dir_path + '\\data\\db_out.json', orient='table', typ='series')
    # clean_data = data_cleaner.clean_from_database(dirty_data['data'])
    clean_data = json.load(open('data\\clean_db_data.json'))
    processing.processing(clean_data, 100, 'b', 1)

    return clean_data


def run():
    start_time = time.time()
    # update_data()
    # count_db_entries()
    get_good_data()
    # process_existing()
    end_time = time.time()
    print("Elapsed time: %g seconds" % (end_time - start_time))


if __name__ == '__main__':
    run()
