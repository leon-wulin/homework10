#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>
#include <chrono>
#include <algorithm>
#include "image.h"
#include "priority_queue.h"

// ===================================================================================================

// distance field method functions
double NaiveDistanceFieldMethod(Image<Color> &input, Image<DistancePixel> &distance_image);
double ImprovedDistanceFieldMethod(Image<Color> &input, Image<DistancePixel> &distance_image);
double FastMarchingMethod(Image<Color> &input, Image<DistancePixel> &distance_image);

// visualization style helper functions
Color Rainbow(double distance, double max_distance);
Color GreyBands(double distance, double max_distance, int num_bands);
Color GradientColor(double distance, double max_distance);

// ===================================================================================================

int main(int argc, char* argv[]) {
  if (argc != 5) {
    std::cerr << "Usage: " << argv[0] << " input.ppm output.ppm distance_field_method visualization_style" << std::endl;
    exit(1);
  }

  // open the input image
  Image<Color> input;
  if (!input.Load(argv[1])) {
    std::cerr << "ERROR: Cannot open input file: " << argv[1] << std::endl;
    exit(1);
  }

  // a place to write the distance values
  Image<DistancePixel> distance_image;
  distance_image.Allocate(input.Width(),input.Height());

  // Start the timer
  std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now(); 

  // calculate the distance field (each function returns the maximum distance value)
  double max_distance = 0;
  if (std::string(argv[3]) == std::string("naive_method")) {
    max_distance = NaiveDistanceFieldMethod(input,distance_image);
  } else if (std::string(argv[3]) == std::string("improved_method")) {
    max_distance = ImprovedDistanceFieldMethod(input,distance_image);
  } else if (std::string(argv[3]) == std::string("pq_with_map")) {
    max_distance = FastMarchingMethod(input,distance_image);
  } else if (std::string(argv[3]) == std::string("pq_with_hash_table")) {
    // EXTRA CREDIT: implement FastMarchingMethod with a hash table
  } else {
    std::cerr << "ERROR: Unknown distance field method: " << argv[3] << std::endl;
    exit(1);
  }
  // Stop the timer
  std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now(); 

  // Calculate the elapsed time in milliseconds
  long long elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count(); 

  std::cout << "The function took " << elapsed_time << " milliseconds to run." << std::endl; 

  // convert distance values to a visualization
  Image<Color> output;
  output.Allocate(input.Width(),input.Height());
  for (int i = 0; i < input.Width(); i++) {
    for (int j = 0; j < input.Height(); j++) {
      double v = distance_image.GetPixel(i,j).getValue();
      if (std::string(argv[4]) == std::string("greyscale")) {
	output.SetPixel(i,j,GreyBands(v,max_distance*1.01,1));
      } else if (std::string(argv[4]) == std::string("grey_bands")) {
	output.SetPixel(i,j,GreyBands(v,max_distance,4));
      } else if (std::string(argv[4]) == std::string("rainbow")) {
	output.SetPixel(i,j,Rainbow(v,max_distance));
      } else if (std::string(argv[4]) == std::string("gradient")) {
    output.SetPixel(i, j, GradientColor(v, max_distance)); 
      } else {
	// EXTRA CREDIT: create other visualizations 
	std::cerr << "ERROR: Unknown visualization style" << std::endl;
	exit(0);
      }
    }
  }
  // save output
  if (!output.Save(argv[2])) {
    std::cerr << "ERROR: Cannot save to output file: " << argv[2] << std::endl;
    exit(1);
  }

  return 0;
}

// ===================================================================================================

double NaiveDistanceFieldMethod(Image<Color> &input, Image<DistancePixel> &distance_image) {
  int w = input.Width();
  int h = input.Height();
  // return the maximum distance value
  double answer = 0;
  // loop over the pixels in the input image
  for (int i = 0; i < w; i++)  {
    for (int j = 0; j < h; j++) {
      double closest = -1;      
      // loop over all other pixels in the input image
      for (int i2 = 0; i2 < w; i2++)  {
	for (int j2 = 0; j2 < h; j2++) {
	  const Color& c = input.GetPixel(i2,j2);      
	  // skip all pixels that are not black
	  if (!c.isBlack()) continue;
	  // calculate the distance between the two pixels
	  double distance = sqrt((i-i2)*(i-i2) + (j-j2)*(j-j2));
	  // store the closest distance to a black pixel
	  if (closest < 0 || distance < closest) {
	    closest = distance;
	  }
	}
      }
      assert (closest >= 0);
      answer = std::max(answer,closest);
      // save the data to the distance image
      DistancePixel& p = distance_image.GetPixel(i,j);
      p.setValue(closest);
    }
  }
  return answer;
}

double ImprovedDistanceFieldMethod(Image<Color> &input, Image<DistancePixel> &distance_image) {

  //
  // IMPLEMENT THIS FUNCTION
  //
  // a small improvement on the NaiveDistanceFieldMethod
  //
    int w = input.Width();
    int h = input.Height();
    double max_distance = 0;

    // Initialize distance values based on black pixels
    for (int i = 0; i < w; i++) {
        for (int j = 0; j < h; j++) {
            const Color& c = input.GetPixel(i, j);
            DistancePixel& p = distance_image.GetPixel(i, j);

            if (c.isBlack()) {
                p.setValue(0);
            }
            else {
                p.setValue(std::numeric_limits<int>::max());
            }
        }
    }

    // Keep updating the distance values until no changes are made
    bool has_changes = true;
    while (has_changes) {
        has_changes = false;
        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {
                const Color& c = input.GetPixel(i, j);

                if (!c.isBlack()) {
                    double current_distance = distance_image.GetPixel(i, j).getValue();
                    double min_distance = current_distance;

                    // Check the 8 neighboring pixels
                    for (int dx = -1; dx <= 1; dx++) {
                        for (int dy = -1; dy <= 1; dy++) {
                            int x = i + dx;
                            int y = j + dy;

                            if (x >= 0 && x < w && y >= 0 && y < h) {
                                double new_distance;
                                if (dx == 0 || dy == 0) {
                                    new_distance = 1;
                                }
                                else {
                                    new_distance = 1.4;
                                }

                                // Update the minimum distance if a better one is found
                                double neighbor_distance = distance_image.GetPixel(x, y).getValue(); 
                                min_distance = std::min(min_distance, neighbor_distance + new_distance); 
                            }
                        }
                    }

                    // Update the current pixel's distance value and set the has_changes flag
                    if (min_distance < current_distance) { 
                        distance_image.GetPixel(i, j).setValue(min_distance); 
                        max_distance = std::max(max_distance, min_distance); 
                        has_changes = true;
                    }
                }
            }
        }
    }

    return max_distance; 
}

double FastMarchingMethod(Image<Color> &input, Image<DistancePixel> &distance_image) {
  
  //
  // IMPLEMENT THIS FUNCTION
  //
  // (using the advancing front method, and a priority queue)
  //
    int w = input.Width();
    int h = input.Height();
    distance_image.Allocate(w, h);
    DistancePixel_PriorityQueue queue;
    double max_distance = 0;

    // Initialize the distance image and populate the priority queue with boundary pixels
    for (int row = 0; row < w; row++) {
        for (int column = 0; column < h; column++) {
            const Color& c = input.GetPixel(row, column);
            DistancePixel& p = distance_image.GetPixel(row, column);

            // Set the correct row and column values for each DistancePixel object
            p.setX(row);
            p.setY(column);

            if (c.isBlack()) {
                p.setValue(0);
                queue.push(&p);
            }
            else {
                p.setValue(std::numeric_limits<double>::max());
            }
        }
    }

    // Fast Marching Method
    while (!queue.empty()) {
        DistancePixel* current = const_cast<DistancePixel*>(queue.top());
        queue.pop();

        // Iterate through neighboring pixels
        for (int row_offset = -1; row_offset <= 1; row_offset++) {
            for (int col_offset = -1; col_offset <= 1; col_offset++) {
                // Skip the current pixel
                if (row_offset == 0 && col_offset == 0) {
                    continue;
                }

                int neighbor_row = current->getX() + row_offset;
                int neighbor_col = current->getY() + col_offset;

                // Check if the neighbor is within the image bounds
                if (neighbor_row >= 0 && neighbor_row < w && neighbor_col >= 0 && neighbor_col < h) {
                    DistancePixel& neighbor = distance_image.GetPixel(neighbor_row, neighbor_col);

                    // Calculate neighbor distance, 1.4 for diagonal neighbors and 1 for others
                    double neighbor_distance = (std::abs(row_offset) == 1 && std::abs(col_offset) == 1) ? 1.4 : 1.0; 

                    // Update the neighbor's distance if a shorter path is found
                    if (current->getValue() + neighbor_distance < neighbor.getValue()) { 
                        neighbor.setValue(current->getValue() + neighbor_distance); 
                        // If the neighbor is already in the priority queue, update its position
                        if (queue.in_heap(&neighbor)) { 
                            queue.update_position(&neighbor); 
                        }
                        // Otherwise, add the neighbor to the priority queue
                        else {
                            queue.push(&neighbor); 
                        }
                    }
                }
            }
        }
    }

    // Find the maximum distance
    for (int row = 0; row < w; row++) {
        for (int column = 0; column < h; column++) {
            const DistancePixel& p = distance_image.GetPixel(row, column);
            if (p.getValue() > max_distance) {
                max_distance = p.getValue();
            }
        }
    }

    return max_distance;
 }

// ===================================================================================================

Color Rainbow(double distance, double max_distance) {
  Color answer;
  if (distance < 0.001) {
    // black
    answer.r = 0; answer.g = 0; answer.b = 0;
  } else if (distance < 0.2*max_distance) {
    // blue -> cyan
    double tmp = distance * 5.0 / max_distance;
    answer.r = 0;
    answer.g = tmp*255;
    answer.b = 255;
  } else if (distance < 0.4*max_distance) {
    // cyan -> green
    double tmp = (distance-0.2*max_distance) * 5.0 / max_distance;
    answer.r = 0;
    answer.g = 255;
    answer.b = (1-tmp*tmp)*255;
  } else if (distance < 0.6*max_distance) {
    // green -> yellow
    double tmp = (distance-0.4*max_distance) * 5.0 / max_distance;
    answer.r = sqrt(tmp)*255;
    answer.g = 255;
    answer.b = 0;
  } else if (distance < 0.8*max_distance) {
    // yellow -> red
    double tmp = (distance-0.6*max_distance) * 5.0 / max_distance;
    answer.r = 255;
    answer.g = (1-tmp*tmp)*255;
    answer.b = 0;
  } else if (distance < max_distance) {
    // red -> white
    double tmp = (distance-0.8*max_distance) * 5.0 / max_distance;
    answer.r = 255;
    answer.g = tmp*255;
    answer.b = tmp*255;
  } else {
    // white
    answer.r = answer.g = answer.b = 255;
  }  
  return answer;
}

Color GreyBands(double distance, double max_value, int num_bands) {
  Color answer;
  if (distance < 0.001) {
    // red
    answer.r = 255; answer.g = 0; answer.b = 0;
  } else {
    // shades of grey
    answer.r = answer.g = answer.b = int(num_bands*256*distance/double(max_value)) % 256;
  }  
  return answer;
}

// This function converts HSV color to RGB color
Color HSVtoRGB(double h, double s, double v) {
    double c = v * s;
    double hh = h / 60.0;
    double x = c * (1 - std::abs(fmod(hh, 2) - 1));
    double m = v - c;

    double r1, g1, b1;

    if (0 <= hh && hh <= 1) {
        r1 = c;
        g1 = x;
        b1 = 0;
    }
    else if (1 <= hh && hh <= 2) {
        r1 = x;
        g1 = c;
        b1 = 0;
    }
    else if (2 <= hh && hh <= 3) {
        r1 = 0;
        g1 = c;
        b1 = x;
    }
    else if (3 <= hh && hh <= 4) {
        r1 = 0;
        g1 = x;
        b1 = c;
    }
    else if (4 <= hh && hh <= 5) {
        r1 = x;
        g1 = 0;
        b1 = c;
    }
    else {
        r1 = c;
        g1 = 0;
        b1 = x;
    }

    Color rgb;
    rgb.r = static_cast<unsigned char>((r1 + m) * 255);
    rgb.g = static_cast<unsigned char>((g1 + m) * 255);
    rgb.b = static_cast<unsigned char>((b1 + m) * 255);

    return rgb;
}

// This function generates a gradient color visualization of the distance field.
// It uses the HSV color space to create a smooth gradient, where the hue value
// is mapped to the normalized distance value.
Color GradientColor(double distance, double max_distance) { 

    if (distance == 0) { 
        Color white; 
        white.r = 255; 
        white.g = 255; 
        white.b = 255; 
        return white; 
    }

    // Normalize the distance value between 0 and 1
    double normalized_distance = distance / max_distance; 

    // Map the normalized distance value to the hue value (0 to 360 degrees)
    double hue = normalized_distance * 360.0; 

    // Set the saturation and value (brightness) to 1 for a full-color gradient
    double saturation = 1.0;
    double value = 1.0;

    // Convert the HSV color to RGB color using the HSVtoRGB function
    Color rgb_color = HSVtoRGB(hue, saturation, value); 

    // Return the gradient color
    return rgb_color; 
}

// ===================================================================================================
