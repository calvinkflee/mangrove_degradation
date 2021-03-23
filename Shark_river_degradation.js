Map.centerObject(ROI, 12);
Map.setOptions('SATELLITE');

// Modelling variables
var nTrees = 50;
var treeSplit = 0; // 0 is the default: sqrt of nPredictors 
var bagFrac = 0.5;
var minLeaf = 1;
var seed = 0;

var GMW2016_image = GMW2016
  .reduceToImage({
    properties: ['pxlval'],
    reducer: ee.Reducer.first()
});

GMW2016_image = GMW2016_image.unmask();
var GMW_mask = ee.Image(1).updateMask(GMW2016_image.neq(1));

// Cloud and shadow masking to null
function cloudMaskingToNull(image){
  var cmask = image.select('pixel_qa').bitwiseAnd(32).eq(0);
  var csmask = image.select('pixel_qa').bitwiseAnd(8).eq(0);
  var masked = image.updateMask(cmask).updateMask(csmask);
  return masked;
}

// Collections
var L8_2015 = L8_SR
  .filterBounds(ROI)
  .filterDate('2015-06-01', '2016-05-31')
  .select(['B2','B3','B4','B5','B6','B7','pixel_qa'], 
    ['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'pixel_qa'])
  .map(cloudMaskingToNull);

Map.addLayer(L8_2015.reduce(ee.Reducer.median()), 
  {bands: ['red_median', 'green_median', 'blue_median'],
    min: 100, max: 800}, '2015 RGB median');

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


var L8_pre = L8_SR
  .filterBounds(ROI)
  .filterDate('2016-06-01', '2017-05-31')
  .select(['B2','B3','B4','B5','B6','B7','pixel_qa'], 
    ['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'pixel_qa'])
  .map(cloudMaskingToNull);

Map.addLayer(L8_pre.reduce(ee.Reducer.median()), 
  {bands: ['red_median', 'green_median', 'blue_median'],
    min: 100, max: 800}, 'pre-Irma RGB median');

var L8_post = L8_SR
  .filterBounds(ROI)
  .filterDate('2017-11-01', '2018-10-31')
  .select(['B2','B3','B4','B5','B6','B7','pixel_qa'], 
    ['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'pixel_qa'])
  .map(cloudMaskingToNull);

Map.addLayer(L8_post.reduce(ee.Reducer.median()), 
  {bands: ['red_median', 'green_median', 'blue_median'],
    min: 100, max: 800}, 'post-Irma RGB median');


var L8_post = L8_post
  .map(cloudMaskingToNull)
  .map(addNDVI)
  .map(addNDWI)
  .map(addLSWI);
  
var post_NDVI = L8_post.select('NDVI');
var post_NDWI = L8_post.select('NDWI');
var post_LSWI = L8_post.select('LSWI');

var avg_NDVI_post = post_NDVI.reduce(ee.Reducer.mean());
var sd_NDVI_post = post_NDVI.reduce(ee.Reducer.stdDev());
var min_NDVI_post = post_NDVI.reduce(ee.Reducer.min());
var max_NDVI_post = post_NDVI.reduce(ee.Reducer.max());
var avg_NDWI_post = post_NDWI.reduce(ee.Reducer.mean());
var sd_NDWI_post = post_NDWI.reduce(ee.Reducer.stdDev());
var avg_LSWI_post = post_LSWI.reduce(ee.Reducer.mean());
var sd_LSWI_post = post_LSWI.reduce(ee.Reducer.stdDev());
var min_LSWI_post = post_LSWI.reduce(ee.Reducer.min());
var max_LSWI_post = post_LSWI.reduce(ee.Reducer.max());


var covariates_post = avg_NDVI_post
  .addBands(sd_NDVI_post)
  .addBands(avg_NDWI_post)
  .addBands(avg_LSWI_post)
  .addBands(sd_LSWI_post)
  .clip(ROI)
  .updateMask(GMW2016_image)

var L8_pre = L8_pre
  .map(cloudMaskingToNull)
  .map(addNDVI)
  .map(addNDWI)
  .map(addLSWI);
  
var pre_NDVI = L8_pre.select('NDVI');
var pre_NDWI = L8_pre.select('NDWI');
var pre_LSWI = L8_pre.select('LSWI');

var avg_NDVI_pre = pre_NDVI.reduce(ee.Reducer.mean());
var sd_NDVI_pre = pre_NDVI.reduce(ee.Reducer.stdDev());
var min_NDVI_pre = pre_NDVI.reduce(ee.Reducer.min());
var max_NDVI_pre = pre_NDVI.reduce(ee.Reducer.max());
var avg_NDWI_pre = pre_NDWI.reduce(ee.Reducer.mean());
var sd_NDWI_pre = pre_NDWI.reduce(ee.Reducer.stdDev());
var avg_LSWI_pre = pre_LSWI.reduce(ee.Reducer.mean());
var sd_LSWI_pre = pre_LSWI.reduce(ee.Reducer.stdDev());
var min_LSWI_pre = pre_LSWI.reduce(ee.Reducer.min());
var max_LSWI_pre = pre_LSWI.reduce(ee.Reducer.max());

var covariates_pre = avg_NDVI_pre
  .addBands(sd_NDVI_pre)
  .addBands(avg_NDWI_pre)
  .addBands(avg_LSWI_pre)
  .addBands(sd_LSWI_pre)
  .clip(ROI)
  .updateMask(GMW2016_image)

// Masking out non-mangrove areas according to GMW2016
Map.addLayer(GMW_mask, {}, 'GMW2016 mask');

Map.addLayer(covariates_pre);
Map.addLayer(covariates_post);

var intact_sample = covariates_pre.sampleRegions({
  collection: intact,
  properties: ['class'],
  scale: 30
});

print(intact_sample);

var degraded_sample = covariates_post.sampleRegions({
  collection: degraded,
  properties: ['class'],
  scale: 30
});

print(degraded_sample);

var training_sample = intact_sample.merge(degraded_sample);
print(training_sample);


// Training models
var RFmodel = ee.Classifier.smileRandomForest(nTrees);

var trainedClassifier = RFmodel.train({
  features: training_sample,
  classProperty: 'class',
  inputProperties: ['NDVI_mean', 
  'NDVI_stdDev',
  'LSWI_mean',
  'LSWI_stdDev',
  'NDWI_mean']
});


// Variable importance code from: https://code.earthengine.google.com/c1e015bc57bec53ba31d2cdaf1e9b40f
var dict = trainedClassifier.explain();

var variable_importance = ee.Feature(null, ee.Dictionary(dict).get('importance'));

var chart =
  ui.Chart.feature.byProperty(variable_importance)
    .setChartType('ColumnChart')
    .setOptions({
      title: 'Random Forest Variable Importance',
      legend: {position: 'none'},
      hAxis: {title: 'Landsat Bands'},
      vAxis: {title: 'Importance'}
    });
print(chart);


// Classify the two pre- and post- images
var plasma = ["0d0887", "3d049b", "6903a5", "8d0fa1", "ae2891", 
  "cb4679", "df6363", "f0844c", "faa638", "fbcc27", "f0f921"];

var output_pre = covariates_pre.classify(trainedClassifier);
Map.addLayer(output_pre, {palette: plasma, min: 0, max: 1},
  'pre_Irma degradation probability');
var output_post = covariates_post.classify(trainedClassifier);
Map.addLayer(output_post, {palette: plasma, min: 0, max: 1},
  'post-Irma degradation');

// Exporting classifications
Export.image.toDrive({
  image: output_pre,
  description: 'pre-Irma_degrad',
  folder: 'Shark_River',
  fileNamePrefix: 'Degradation_binary_pre_Irma',
  scale: 30,
  region: ROI,
  crs: 'EPSG:32617'
});

Export.image.toDrive({
  image: output_post,
  description: 'post-Irma_degrad',
  folder: 'Shark_River',
  fileNamePrefix: 'Degradation_binary_post_Irma',
  scale: 30,
  region: ROI,
  crs: 'EPSG:32617'
});

