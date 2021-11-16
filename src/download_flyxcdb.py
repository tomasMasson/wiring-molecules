#!/usr/bin/env python3

"Download Fly Extracellular Database (FlyXCDB) table, published in the Journal of Molecular Biology"

import requests
import pandas as pd
from bs4 import BeautifulSoup

# Define the location of the data
URL = 'http://prodata.swmed.edu/FlyXCDB/info.list.new21_26.html'

# Retrieve the data using a request
r = requests.get(URL)

# Get the content in HTML format
html = r.text

# Parse the HTML file with BS4
soup = BeautifulSoup(html, 'html.parser')

# Search for any kind of tabular data (ignores descriptions and supporting data)
table = soup.find_all('table')

# Define the table that we are interested (in this case, we have only one table in the HTML file)
table = table[0]

# Extract the columns headers for the CSV table
headers = [th.text.strip()
           for th in table.find('tr').find_all('th')]

# Now, extract the data fields for each row
rows = []
# First, iterate over all the record rows
for tr in table.find_all('tr')[1:]:
    # Then, create a list to put all the individual data fields
    cells = []
    # Find these data fields
    tds = tr.find_all('td')
    # Iterate over the data fields and append them to cells (remove trailing characters)
    for td in tds:
        cells.append(td.text.strip())
    # Finally, add all the extracted data to the rows list
    rows.append(cells)

# Create a Pandas DataFrame using the columns and rows list
df = pd.DataFrame(rows, columns=headers)
# Save the DataFrame to a CSV file
df.to_csv('flyxcdb.csv', index=False)
