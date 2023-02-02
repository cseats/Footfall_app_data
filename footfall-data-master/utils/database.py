import os
import pyodbc
import pandas as pd
import time
import plotly.express as px

# from utils import graph_data

dir_path = os.path.dirname(os.path.realpath(__file__))

def graph_data(df):
    cnt = 0
    datecnt = 0
    date_data = dict()
    dates_ = []
    occurances = []

    clean_date = []
    user = []
    print(df.columns)
    for i in range(len(df["CreatedAt"])):
        if df["CreatedAt"][i][0:4] == "2022":
            print(df["CreatedAt"][i][0:4])
            clean_date.append(df["CreatedAt"][i][0:10])
            user.append(df["DeviceId"][i])

    df_clean = pd.DataFrame(list(zip(clean_date,user)),columns=["Date","User"])

    #determine Uses by date
    uniq_dates = df_clean.Date.unique()
    print(uniq_dates)
    date_oc = dict()
    for i in uniq_dates:
        num_date = df_clean.Date.str.contains(i).sum()
        date_oc[i] = num_date

    print(date_oc)
    print(date_oc.keys())
    print(date_oc.values())
    date_df = pd.DataFrame(list(zip(list(date_oc.keys()),list(date_oc.values()))),columns=["Date","Uses"])
    fig = px.bar(date_df, x='Date', y='Uses')
    fig.show()

    uniq_user = df_clean.User.unique()
    print(df_clean)
    user_use = dict()
    for i in uniq_user:
        user_num = df_clean.User.str.contains(i).sum()
        user_use[i] = user_num
    print(user_use)
    usere_df = pd.DataFrame(list(zip(list(user_use.keys()),list(user_use.values()))),columns=["Users","Uses"])
    usere_df = usere_df.sort_values(by=["Uses"])
    print(usere_df)
    fig = px.bar(usere_df, x='Users', y='Uses')
    fig.show()
    print(f"There are {len(uniq_user)} users durring this time period")


    #determine Uses by date
#
#         cur_date = i[0:10]
#         print(cur_date)
#         if cnt ==0:
#             ind_date = cur_date
#
#         if cur_date==ind_date:
#             datecnt = datecnt+1
#             date_data[ind_date] = datecnt
#
#         else:
#             dates_.append(ind_date)
#             occurances.append(datecnt)
#             datecnt = 0
#             cnt = 0
#             ind_date = cur_date
#         cnt +=1
#
#     print(date_data)
#     date_df = pd.DataFrame.from_dict(date_data, orient='index',columns=["occur"])
#     print(list(date_df.index))
#     print(date_df["occur"].tolist())
#
#
#
# # pd.DataFrame(list(zip(lst, lst2))
#     date_df2 = pd.DataFrame(list(zip(list(date_df.index), date_df["occur"].tolist())),columns=["Date","Occurance"])
#     print(date_df)
#     print(date_df2)
#     # import plotly.express as px
#     fig = px.bar(date_df2, x='Date', y='Occurance')
#     fig.show()
#     return i


def connect():
    conn = pyodbc.connect(
        r'DRIVER={SQL Server};'
        r'SERVER=footfall.database.windows.net;'
        r'DATABASE=arupfootfallwa209;'
        r'UID=dbreader;'
        r'PWD=re@der2O!7'
                # r'DRIVER={SQL Server};'
                # r'SERVER=footfall.database.windows.net;'
                # r'DATABASE=arupfootfallwa209;'
                # r'UID=dbreader;'
                # r'PWD=re@der2O!7'
    )
    return conn


def close(conn):
    conn.close()


def run_pandas_query(conn, query):
    df = pd.read_sql(query, conn)
    close(conn)
    return df


def run_query(conn, query):
    cursor = conn.cursor()
    res = cursor.execute(query).fetchall()[0][0]
    close(conn)
    return res


if __name__ == '__main__':
    start_time = time.time()

    connection = connect()
    get_all_query = """
      SELECT TOP (2000) [Id]
          ,RecordedData
          ,DeviceId
          ,CreatedAt
          ,Deleted
          ,MeasuredPeakResponseFactor
          ,MeasuredNaturalFreq
          ,Metadata
          ,Location
      FROM dbo.VibrationItems
      order by [CreatedAt] desc
                   """
    df = run_pandas_query(connection,get_all_query)
    print(df)
    print(df["CreatedAt"])

    graph_data(df)


    # res = run_query(connection, get_all_query)

    # connection.commit()
    # print(res)
    # print(get_all_query)
    # close(connection)
    end_time = time.time()
    print("Elapsed time was %g seconds" % (end_time - start_time))
