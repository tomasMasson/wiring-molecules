#!/usr/bin/env python3

"Download Drosophila melanogaster extracellular domain batabase (FlyXCDB) table, published in the Journal of Molecular Biology"

import click
import requests
import pandas as pd
from bs4 import BeautifulSoup

def scrape_url(url):
    """
    Scrap the content from the input URL using Beautiful Soup.
    """

    # Retrieve the data using a request
    r = requests.get(url)

    # Get the content in HTML format
    html = r.text

    # Parse the HTML file with BS4
    soup = BeautifulSoup(html, 'html.parser')
    
    # Print a brief report
    print('FlyXCDB was succesfully scraped')

    return soup

def extract_flyxcdb_data(html, output):
    """
    Parse the HTML data from FlyXCDB and write the table into a the output file.
    """

    # Search for any kind of tabular data (ignores descriptions and supporting data)
    table = html.find_all('table')

    # Define the table that we are interested (in this case, we have only one table)
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
    df.to_csv(output, index=False)

    print(f"FlyXCDB data has been saved to '{output}'")

# FlyXCDB data URL
URL = 'http://prodata.swmed.edu/FlyXCDB/info.list.new21_26.html'

# CLI options
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--html',
              default=URL,
              help="Url used to scrap the data from FlyXCDB. Default is FlyXCDB URL")
@click.option('--output',
              default='flyxcdb_data.csv',
              help="Output name for the scraped table. Defaults is 'flyxcdb_data.csv'")

# CLI main function
def command_line_interface(html, output):
    """
    Provides the HTML address and the output file name through the CLI
    """

    content = scrape_url(html)

    extract_flyxcdb_data(content, output)

if __name__ == '__main__':
    command_line_interface()
