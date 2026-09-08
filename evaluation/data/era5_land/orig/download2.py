import cdsapi

dataset = "reanalysis-era5-land"
request = {
    "variable": [
        "snow_density",
        "snow_depth_water_equivalent"
    ],
    "year": "2025",
    "month": [
        "09", "10", "11",
        "12"
    ],
    "day": [
        "01", "02", "03",
        "04", "05", "06",
        "07", "08", "09",
        "10", "11", "12",
        "13", "14", "15",
        "16", "17", "18",
        "19", "20", "21",
        "22", "23", "24",
        "25", "26", "27",
        "28", "29", "30",
        "31"
    ],
    "time": ["00:00"],
    "data_format": "netcdf",
    "download_format": "unarchived"
}

filename="ERA5land_00Z_1Sep_31Dec25.nc"
client = cdsapi.Client()
client.retrieve(dataset, request, filename)
