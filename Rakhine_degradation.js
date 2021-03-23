Map.setOptions('SATELLITE');

site = rakhineMangroves;

Map.centerObject(site);

var outFolder = 'Rakhine';
var nTrees = 50;


var GMW2016_image = GMW2016
  .reduceToImage({
    properties: ['pxlval'],
    reducer: ee.Reducer.first()
});

GMW2016_image = GMW2016_image.unmask();
var GMW_mask = ee.Image(1).updateMask(GMW2016_image.neq(1));

function createYearBand(image) {
    // add a new band with year to all images
    var year = ee.Date(image.get('system:time_start')).get('year').subtract(1984);
    return ee.Image(year).byte().addBands(image);
}

var createTimeBand = function(image) {
  // Scale milliseconds by a large constant to avoid very small slopes
  // in the linear regression output. 3e11 scales VERY roughly to years.
  return image.addBands(image.metadata('system:time_start').divide(3e11));
};

// Cloud and shadow masking to null
function cloudMaskingToNull(image){
  var cmask = image.select('pixel_qa').bitwiseAnd(32).eq(0);
  var csmask = image.select('pixel_qa').bitwiseAnd(8).eq(0);
  var masked = image.updateMask(cmask).updateMask(csmask);
  return masked;
}

// calculate NDVI and add to collection
function addNDVI(image) {
  return image.addBands(image.normalizedDifference(['nir', 'red']).rename('NDVI'));
}


// calculate LSWI and add to collection
function addLSWI(image){
  return image.addBands(image.normalizedDifference(['nir', 'swir1']).rename('LSWI'));
}

// NDMI are the same as LSWI (there's also a NDWI (Gao 1996) that's the same).

// Calculate NDWI and add to collection
function addNDWI(image){
  return image.addBands(image.normalizedDifference(['green', 'nir']).rename('NDWI'));
}

var L8_2016 = L8_SR
  .filterBounds(site)
  .filterDate('2016-01-01', '2016-12-31')
  .select(['B2','B3','B4','B5','B6','B7','pixel_qa'], 
    ['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'pixel_qa']);

Map.addLayer(L8_2016.reduce(ee.Reducer.median()), 
  {bands: ['red_median', 'green_median', 'blue_median'],
    min: 0, max: 2500}, '2016 RGB median');

var L8_col = L8_SR
  .filterBounds(site)
  .filterDate('2017-01-01', '2017-12-31')
  .select(['B2','B3','B4','B5','B6','B7','pixel_qa'], 
    ['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'pixel_qa']);
    

L8_col = L8_col
  .map(cloudMaskingToNull)
  .map(addNDVI)
  .map(addNDWI)
  .map(addLSWI);
  
var L8_NDVI = L8_col.select('NDVI');
var L8_NDWI = L8_col.select('NDWI');
var L8_LSWI = L8_col.select('LSWI');

var avg_NDVI = L8_NDVI.reduce(ee.Reducer.mean());
var sd_NDVI = L8_NDVI.reduce(ee.Reducer.stdDev());
var avg_NDWI = L8_NDWI.reduce(ee.Reducer.mean());
var sd_NDWI = L8_NDWI.reduce(ee.Reducer.stdDev());
var avg_LSWI = L8_LSWI.reduce(ee.Reducer.mean());
var sd_LSWI = L8_LSWI.reduce(ee.Reducer.stdDev());

var NDVI_25p = L8_NDVI.reduce(ee.Reducer.percentile([25]));
var NDVI_75p = L8_NDVI.reduce(ee.Reducer.percentile([75]));
var NDVI_IQR = NDVI_75p.subtract(NDVI_25p);


// ALOS PALSAR
function PALSARConversion(image){
  return image.pow(2).log10().multiply(ee.Image(10)).subtract(ee.Image(83.0));
}

function refinedLee(img) {
  // img must be in natural units, i.e. not in dB!
  // Set up 3x3 kernels 

  // convert to natural.. do not apply function on dB!
  var myimg = toNatural(img);

  var weights3 = ee.List.repeat(ee.List.repeat(1,3),3);
  var kernel3 = ee.Kernel.fixed(3,3, weights3, 1, 1, false);

  var mean3 = myimg.reduceNeighborhood(ee.Reducer.mean(), kernel3);
  var variance3 = myimg.reduceNeighborhood(ee.Reducer.variance(), kernel3);

  // Use a sample of the 3x3 windows inside a 7x7 windows to determine gradients and directions
  var sample_weights = ee.List([[0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0], [0,1,0,1,0,1,0], [0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0]]);

  var sample_kernel = ee.Kernel.fixed(7,7, sample_weights, 3,3, false);

  // Calculate mean and variance for the sampled windows and store as 9 bands
  var sample_mean = mean3.neighborhoodToBands(sample_kernel); 
  var sample_var = variance3.neighborhoodToBands(sample_kernel);

  // Determine the 4 gradients for the sampled windows
  var gradients = sample_mean.select(1).subtract(sample_mean.select(7)).abs();
  gradients = gradients.addBands(sample_mean.select(6).subtract(sample_mean.select(2)).abs());
  gradients = gradients.addBands(sample_mean.select(3).subtract(sample_mean.select(5)).abs());
  gradients = gradients.addBands(sample_mean.select(0).subtract(sample_mean.select(8)).abs());

  // And find the maximum gradient amongst gradient bands
  var max_gradient = gradients.reduce(ee.Reducer.max());

  // Create a mask for band pixels that are the maximum gradient
  var gradmask = gradients.eq(max_gradient);

  // duplicate gradmask bands: each gradient represents 2 directions
  gradmask = gradmask.addBands(gradmask);

  // Determine the 8 directions
  var directions = sample_mean.select(1).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(7))).multiply(1);
  directions = directions.addBands(sample_mean.select(6).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(2))).multiply(2));
  directions = directions.addBands(sample_mean.select(3).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(5))).multiply(3));
  directions = directions.addBands(sample_mean.select(0).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(8))).multiply(4));
  // The next 4 are the not() of the previous 4
  directions = directions.addBands(directions.select(0).not().multiply(5));
  directions = directions.addBands(directions.select(1).not().multiply(6));
  directions = directions.addBands(directions.select(2).not().multiply(7));
  directions = directions.addBands(directions.select(3).not().multiply(8));

  // Mask all values that are not 1-8
  directions = directions.updateMask(gradmask);

  // "collapse" the stack into a singe band image (due to masking, each pixel has just one value (1-8) in it's directional band, and is otherwise masked)
  directions = directions.reduce(ee.Reducer.sum());  

  var sample_stats = sample_var.divide(sample_mean.multiply(sample_mean));

  // Calculate localNoiseVariance
  var sigmaV = sample_stats.toArray().arraySort().arraySlice(0,0,5).arrayReduce(ee.Reducer.mean(), [0]);

  // Set up the 7*7 kernels for directional statistics
  var rect_weights = ee.List.repeat(ee.List.repeat(0,7),3).cat(ee.List.repeat(ee.List.repeat(1,7),4));

  var diag_weights = ee.List([[1,0,0,0,0,0,0], [1,1,0,0,0,0,0], [1,1,1,0,0,0,0], 
    [1,1,1,1,0,0,0], [1,1,1,1,1,0,0], [1,1,1,1,1,1,0], [1,1,1,1,1,1,1]]);

  var rect_kernel = ee.Kernel.fixed(7,7, rect_weights, 3, 3, false);
  var diag_kernel = ee.Kernel.fixed(7,7, diag_weights, 3, 3, false);

  // Create stacks for mean and variance using the original kernels. Mask with relevant direction.
  var dir_mean = myimg.reduceNeighborhood(ee.Reducer.mean(), rect_kernel).updateMask(directions.eq(1));
  var dir_var = myimg.reduceNeighborhood(ee.Reducer.variance(), rect_kernel).updateMask(directions.eq(1));

  dir_mean = dir_mean.addBands(myimg.reduceNeighborhood(ee.Reducer.mean(), diag_kernel).updateMask(directions.eq(2)));
  dir_var = dir_var.addBands(myimg.reduceNeighborhood(ee.Reducer.variance(), diag_kernel).updateMask(directions.eq(2)));

  // and add the bands for rotated kernels
  for (var i=1; i<4; i++) {
    dir_mean = dir_mean.addBands(myimg.reduceNeighborhood(ee.Reducer.mean(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
    dir_var = dir_var.addBands(myimg.reduceNeighborhood(ee.Reducer.variance(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
    dir_mean = dir_mean.addBands(myimg.reduceNeighborhood(ee.Reducer.mean(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
    dir_var = dir_var.addBands(myimg.reduceNeighborhood(ee.Reducer.variance(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
  }

  // "collapse" the stack into a single band image (due to masking, each pixel has just one value in it's directional band, and is otherwise masked)
  dir_mean = dir_mean.reduce(ee.Reducer.sum());
  dir_var = dir_var.reduce(ee.Reducer.sum());

  // A finally generate the filtered value
  var varX = dir_var.subtract(dir_mean.multiply(dir_mean).multiply(sigmaV)).divide(sigmaV.add(1.0));

  var b = varX.divide(dir_var);

  var result = dir_mean.add(b.multiply(myimg.subtract(dir_mean)));
  // return(result);
  return(ee.Image(toDB(result.arrayGet(0))).rename("filter"));
 
}



var PALSAR_2017 = PALSAR
  .filterDate('2017-01-01', '2017-12-31');
  
var PALSAR_2017_HH = PALSARConversion(PALSAR_2017.select('HH').first());
var PALSAR_2017_HV = PALSARConversion(PALSAR_2017.select('HV').first());
var PALSAR_2017_ratio = PALSAR_2017_HH.divide(PALSAR_2017_HV);

var PALSAR_2017_RGB = PALSAR_2017_HH
  .addBands(PALSAR_2017_HV)
  .addBands(PALSAR_2017_ratio);

var PALSAR_2017_HH_filtered = refinedLee(PALSAR_2017_HH).rename('HH');
var PALSAR_2017_HV_filtered = refinedLee(PALSAR_2017_HV).rename('HV');


// Masking out non-mangrove areas according to GMW2016
Map.addLayer(GMW_mask, {}, 'GMW2016 mask');

// Selecting variables.
var covariates = avg_NDVI
  .addBands(PALSAR_2017_HH_filtered)
  .addBands(PALSAR_2017_HV_filtered)
  .addBands(sd_NDVI)
  .addBands(NDVI_IQR)
  .addBands(avg_NDWI)
  .addBands(avg_LSWI)
  .addBands(sd_LSWI)
  .clip(site)
  .updateMask(GMW2016_image);

// Three class sequential training
// Intact vs not intact (will be similar to current model)
var intact_points = ee.FeatureCollection([
  ee.Feature(degraded, {'class': 0}),
  ee.Feature(collapsed, {'class': 0}),
  ee.Feature(intact, {'class': 1}),
]);

print(intact_points);

// Creating and training sample (actual values)
var intact_sample = covariates.sampleRegions({
  collection: intact_points,
  properties: ['class'],
  scale: 30
});
print('intact sample', intact_sample);

// Splitting training and testing datasets
var intact_withRandom = intact_sample.randomColumn('random');

var split = 0.8;  // Roughly 80% training, 20% testing.
var intact_trainingPartition = intact_withRandom.filter(ee.Filter.lt('random', split));
var intact_testingPartition = intact_withRandom.filter(ee.Filter.gte('random', split));

// Training models
var RFmodel = ee.Classifier.smileRandomForest(nTrees);

// Trained with 80% of our data.
var intact_trainedClassifier = RFmodel.train({
  features: intact_trainingPartition,
  classProperty: 'class',
  inputProperties: ['NDVI_mean', 
  'HH',
  'HV',
  'NDVI_stdDev',
  'LSWI_mean',
  'LSWI_stdDev',
  'NDWI_mean']
});

// Variable importance code from: https://code.earthengine.google.com/c1e015bc57bec53ba31d2cdaf1e9b40f
var intact_dict = intact_trainedClassifier.explain();
print('Explain:', intact_dict);

var intact_variable_importance = ee.Feature(null, ee.Dictionary(intact_dict).get('importance'));

var intact_chart =
  ui.Chart.feature.byProperty(intact_variable_importance)
    .setChartType('ColumnChart')
    .setOptions({
      title: 'Random Forest Variable Importance for Intact model',
      legend: {position: 'none'},
      hAxis: {title: 'Covariates'},
      vAxis: {title: 'Importance'}
    });
print(intact_chart);



// Classify the test FeatureCollection.
var intact_test = intact_testingPartition.classify(intact_trainedClassifier);

// Print the confusion matrix.
var intact_confusionMatrix = intact_test.errorMatrix('class', 'classification');
print('Confusion Matrix for Intact model', intact_confusionMatrix);
print('Validation overall accuracy for intact model: ', intact_confusionMatrix.accuracy());


// Classify the whole image
var intact_output = covariates.classify(intact_trainedClassifier);

var plasma = ["0d0887", "3d049b", "6903a5", "8d0fa1", "ae2891", 
  "cb4679", "df6363", "f0844c", "faa638", "fbcc27", "f0f921"];

Map.addLayer(intact_output, {palette: plasma, min: 0, max: 1}, 'Intact mangrove map');




// Collapsed vs not collapsed (I'm expecting this to be very bad.)
var collapsed_points = ee.FeatureCollection([
  ee.Feature(degraded, {'class': 1}),
  ee.Feature(collapsed, {'class': 0}),
  ee.Feature(intact, {'class': 1}),
]);

print(collapsed_points);

// Creating and training sample (actual values)
var collapsed_sample = covariates.sampleRegions({
  collection: collapsed_points,
  properties: ['class'],
  scale: 30
});
print('collapsed sample', collapsed_sample);

// Splitting training and testing datasets
var collapsed_withRandom = collapsed_sample.randomColumn('random');

var split = 0.8;  // Roughly 80% training, 20% testing.
var collapsed_trainingPartition = collapsed_withRandom.filter(ee.Filter.lt('random', split));
var collapsed_testingPartition = collapsed_withRandom.filter(ee.Filter.gte('random', split));

// Training models
var RFmodel = ee.Classifier.smileRandomForest(nTrees);

// Trained with 80% of our data.
var collapsed_trainedClassifier = RFmodel.train({
  features: collapsed_trainingPartition,
  classProperty: 'class',
  inputProperties: ['NDVI_mean', 
  'HH',
  'HV',
  'NDVI_stdDev',
  'LSWI_mean',
  'LSWI_stdDev',
  'NDWI_mean']
});

// Variable importance code from: https://code.earthengine.google.com/c1e015bc57bec53ba31d2cdaf1e9b40f
var collapsed_dict = collapsed_trainedClassifier.explain();
print('Explain:', collapsed_dict);

var collapsed_variable_importance = ee.Feature(null, ee.Dictionary(collapsed_dict).get('importance'));

var collapsed_chart =
  ui.Chart.feature.byProperty(collapsed_variable_importance)
    .setChartType('ColumnChart')
    .setOptions({
      title: 'Random Forest Variable Importance for collapsed model',
      legend: {position: 'none'},
      hAxis: {title: 'Covariates'},
      vAxis: {title: 'Importance'}
    });
print(collapsed_chart);



// Classify the test FeatureCollection.
var collapsed_test = collapsed_testingPartition.classify(collapsed_trainedClassifier);

// Print the confusion matrix.
var collapsed_confusionMatrix = collapsed_test.errorMatrix('class', 'classification');
print('Confusion Matrix for collapsed model', collapsed_confusionMatrix);
print('Validation overall accuracy for collapsed model: ', collapsed_confusionMatrix.accuracy());


// Classify the whole image
var collapsed_output = covariates.classify(collapsed_trainedClassifier);

Map.addLayer(collapsed_output, {palette: plasma, min: 0, max: 1},
  'collapsed mangrove map');

// Combining the two classifications
var three_class = intact_output.add(collapsed_output);

Map.addLayer(three_class, {palette: plasma, min: 0, max: 2},
  'three class mangrove map');
  
Export.image.toDrive({
  image: intact_output,
  description: 'exportIntact',
  folder: outFolder,
  fileNamePrefix: 'intact_binary',
  scale: 30,
  region: site,
  maxPixels: 5e9
});

Export.image.toDrive({
  image: collapsed_output,
  description: 'exportCollapsed',
  folder: outFolder,
  fileNamePrefix: 'collapsed_binary',
  scale: 30,
  region: site,
  maxPixels: 5e9
});


var covariates_output = covariates.toDouble();
print(covariates_output);

Export.image.toDrive({
  image: covariates_output,
  description: 'exportCovariates',
  folder: outFolder,
  fileNamePrefix: 'covariates_binary',
  scale: 30,
  region: site,
  maxPixels: 5e9
});
