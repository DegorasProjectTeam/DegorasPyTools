import json

class StarDataManager:

    def __init__(self, json_file_path):
        self.json_file_path = json_file_path
        self.stars_data = {}

    # Read stars data from the JSON file
    def read_json(self):
        with open(self.json_file_path, 'r') as file:
            self.stars_data = json.load(file)
        return self.stars_data

    # Write stars data to the JSON file
    def write_json(self, filtered_stars, output_json_path):
        self.stars_data['stars'] = filtered_stars
        with open(output_json_path, 'w') as output_file:
            json.dump(self.stars_data, output_file, indent=4)
