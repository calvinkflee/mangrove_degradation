Map.setOptions('SATELLITE');

site = rakhineMangroves;

Map.centerObject(site);

var globFunctions = require('users/calvinkflee/default:calvinFunctions');
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
  .map(globFunctions.cloudMaskingToNull)
  .map(globFunctions.addNDVI)
  .map(globFunctions.addNDWI)
  .map(globFunctions.addSVVI)
  .map(globFunctions.addLSWI);

var L8_NDVI = L8_col.select('NDVI');
var L8_NDWI = L8_col.select('NDWI');
var L8_SVVI = L8_col.select('SVVI');
var L8_LSWI = L8_col.select('LSWI');

var avg_NDVI = L8_NDVI.reduce(ee.Reducer.mean());
var sd_NDVI = L8_NDVI.reduce(ee.Reducer.stdDev());
var avg_NDWI = L8_NDWI.reduce(ee.Reducer.mean());
var sd_NDWI = L8_NDWI.reduce(ee.Reducer.stdDev());
var avg_LSWI = L8_LSWI.reduce(ee.Reducer.mean());
var sd_LSWI = L8_LSWI.reduce(ee.Reducer.stdDev());

var max_SVVI = L8_SVVI.reduce(ee.Reducer.max());

var NDVI_25p = L8_NDVI.reduce(ee.Reducer.percentile([25]));
var NDVI_75p = L8_NDVI.reduce(ee.Reducer.percentile([75]));
var NDVI_IQR = NDVI_75p.subtract(NDVI_25p);


// ALOS PALSAR
var PALSAR_2017 = PALSAR
  .filterDate('2017-01-01', '2017-12-31');

var PALSAR_2017_HH = globFunctions.PALSARConversion(PALSAR_2017.select('HH').first());
var PALSAR_2017_HV = globFunctions.PALSARConversion(PALSAR_2017.select('HV').first());
var PALSAR_2017_ratio = PALSAR_2017_HH.divide(PALSAR_2017_HV);

var PALSAR_2017_RGB = PALSAR_2017_HH
  .addBands(PALSAR_2017_HV)
  .addBands(PALSAR_2017_ratio);

var PALSAR_2017_HH_filtered = globFunctions.refinedLee(PALSAR_2017_HH).rename('HH');
var PALSAR_2017_HV_filtered = globFunctions.refinedLee(PALSAR_2017_HV).rename('HV');


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

Map.addLayer(intact_output, {palette: globFunctions.plasma, min: 0, max: 1},
  'Intact mangrove map');




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

Map.addLayer(collapsed_output, {palette: globFunctions.plasma, min: 0, max: 1},
  'collapsed mangrove map');

// Combining the two classifications
var three_class = intact_output.add(collapsed_output);

Map.addLayer(three_class, {palette: globFunctions.plasma, min: 0, max: 2},
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
