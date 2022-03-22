import requests
import sys
import os

os.mkdir("Data")


names = ["Europe_2019_hot.grib", "Europe_2019_cold_2.grib", "Europe_2020_cold_1.grib", "Europe_2019_hot.txt", "Europe_2019_hot.txt"]

for name in names:

    print("Downloading %s... " % name, end = "")

    url = 'https://ifisc.uib-csic.es/users/alex/data_pd_risk_example/%s' % name

    r = requests.get(url, allow_redirects=True)

    open('Data/%s' % name, 'wb').write(r.content)

    print("Done!")
