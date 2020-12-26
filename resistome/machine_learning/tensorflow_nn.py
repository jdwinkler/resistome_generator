from resistome import constants
import numpy
import scipy.stats
import os
from resistome.graphics import visualization as glib
from resistome.sql import sql_interface
from collections import defaultdict

try:
    import tensorflow as tf
    from sklearn.model_selection import train_test_split
    from sklearn.utils import resample
except ImportError:
    raise


def generate_tf_features_labels(feature_type, category_exclusions):
    """

    Builds a daataset from the Resistome using the specified feature types, excluding any categories not of interest.

    Essentially collates all the feature types for the type requested (all genes mutated in the database, all GO
    processes, etc), converts each genotype into a binary vector representation where 1 = feature mutated, and converts
    all phenotype labels into another binary vector representing 1 = possess a general, category-level phenotype.

    Categories in category_exclusions are ignored when processing.

    Returns set of all features, matrix of binary vectors representing genotypes, dict int => list where int is a 
    category index, and int_to_category_dict is a dict (int : str) to convert back to the category string.


    :param feature_type: 
    :param category_exclusions: 
    :return: 
    """

    from resistome.machine_learning import recommendation_system

    vectors, feature_set, general_categories, category_dict = recommendation_system.generate_vector_set(feature_type)

    feature_set = [x.name for x in feature_set]

    general_categories = general_categories - category_exclusions
    general_categories = sorted(general_categories)

    category_to_int_dict = dict()
    int_to_category_dict = dict()
    for i, category in zip(range(0, len(general_categories)), general_categories):
        category_to_int_dict[category] = i
        int_to_category_dict[i] = category

    # convert into integer vectors
    features_to_int_dict = dict()
    feature_set = sorted(feature_set)
    for i, feature in zip(range(0, len(feature_set)), feature_set):
        features_to_int_dict[feature] = i

    label_data = defaultdict(list)

    indices = set(range(0, len(general_categories)))

    feature_dataset = []

    for feature_vector in vectors:

        vid = feature_vector.id
        integer_features = feature_vector.integer_vector(features_to_int_dict, len(feature_set))

        # 2 here accounts for negative categories-sensitive
        category_int_features = [0] * len(general_categories)
        hit_categories = set()

        for (cat, t) in category_dict[vid]:

            if cat not in general_categories:
                continue

            if category_to_int_dict[cat] in hit_categories:
                continue

            category_int_features[category_to_int_dict[cat]] = 1
            hit_categories.add(category_to_int_dict[cat])
            label_data[category_to_int_dict[cat]].append(1)

        indices_to_update = indices - hit_categories

        for i in indices_to_update:
            label_data[i].append(0)

        feature_dataset.append(integer_features)

    return feature_set, feature_dataset, label_data, int_to_category_dict


def get_tensorflow_test_datasets(feature_set, feature_dataset, label_data, test_size=None, train_size=None):
    """

    Builds input test and train datasets for processing in Tensorflow.

    :param feature_set: 
    :param feature_dataset: 
    :param label_data: 
    :param test_size: 
    :param train_size: 
    :return: 
    """

    feature_columns = [tf.contrib.layers.real_valued_column('test',
                                                            dimension=len(feature_set),
                                                            dtype=tf.int32)]
    test_datasets = dict()
    train_datasets = dict()

    if test_size is None and train_size is None:
        # default split
        test_size = 0.25
        train_size = 0.75
    elif test_size is not None:
        # complement of test size
        train_size = 1.0 - test_size
    else:
        # complement of train size
        test_size = 1.0 - train_size

    for i in label_data:

        X = feature_dataset
        y = label_data[i]

        X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                            test_size=test_size,
                                                            train_size=train_size,
                                                            shuffle=True)
        pos_X = []
        pos_Y = []

        neg_X = []
        neg_Y = []

        for x, y in zip(X_train, y_train):
            if y == 1:
                pos_X.append(x)
                pos_Y.append(y)
            else:
                neg_X.append(x)
                neg_Y.append(y)

        test_datasets[i] = (X_test, y_test)
        # split into pos/neg examples to class re-weighting
        # but also ensure that these are drawn from a different subset than the test examples
        train_datasets[i] = (pos_X, pos_Y, neg_X, neg_Y)

    return feature_columns, test_datasets, train_datasets


def train_tensorflow_model(feature_type, scratch_dir, category_exclusions=set()):
    """

    Trains a set of binary classification DNNs from resistome data, excluding any categories specified in 
    category_exclusions. These DNNs have only a single layer of 32 nodes and are trained with 75% of the Resistome
    dataset and tested with random subsets of the other 25%

    Training data is subsampled such that the number of positive (in-class) and negative (out of class) samples are 
    the same. However, some categories have very few positive examples, such as high_ph or high_gravity, so you
    may wish to exclude them from training.

    An important caveat to keep in mind when judging the accuracy of the classifiers is that phenotypes are not often
    screened in-depth: so it is likely the genotype-phenotype relationships are incomplete for each genotype in the
    database.

    Models are saved under the inputs/neural_network_models folder.

    :param feature_type: 
    :param scratch_dir: 
    :param category_exclusions: 
    :return: 
    """

    if feature_type not in {'go', 'gene', 'metabolite'}:
        raise AssertionError('Unknown feature type: %s' % feature_type)

    feature_set, feature_dataset, label_data, label_names = generate_tf_features_labels(feature_type,
                                                                                        category_exclusions)

    feature_columns, test_datasets, train_datasets = get_tensorflow_test_datasets(feature_set,
                                                                                  feature_dataset,
                                                                                  label_data,
                                                                                  test_size=0.25)
    for i in label_data:

        def train_input():

            (pos_X, pos_Y, neg_X, neg_Y) = train_datasets[i]

            # resample with replacement in case len(pos_X or neg_X)  < 500
            pX, pY = resample(pos_X, pos_Y, n_samples=500)
            nX, nY = resample(neg_X, neg_Y, n_samples=500)

            feature_data = tf.concat([tf.convert_to_tensor(pX, dtype=tf.int32),
                                      tf.convert_to_tensor(nX, dtype=tf.int32)], axis=0)

            label_data = tf.concat([tf.convert_to_tensor(pY, dtype=tf.int32),
                                    tf.convert_to_tensor(nY, dtype=tf.int32)], axis=0)

            feature_dict = {'test': feature_data}

            return feature_dict, label_data

        def test_input():

            (X_test, y_test) = test_datasets[i]

            # random subsample of test data
            _, X_test, _, y_test = train_test_split(X_test,
                                                    y_test,
                                                    test_size=0.25,
                                                    shuffle=True)

            feature_data = tf.convert_to_tensor(X_test, dtype=tf.int32)
            label_data = tf.convert_to_tensor(y_test, dtype=tf.int32)
            feature_dict = {'test': feature_data}

            return feature_dict, label_data

        category_name = label_names[i]

        if category_name == 'quaternary_ammonia':
            continue

        model_dir = os.path.join(constants.INPUT_DIR, 'neural_network_models', category_name)
        new_scratch = os.path.join(scratch_dir, category_name)

        # todo: will break if directories are not empty
        # todo: portability?
        try:
            os.removedirs(new_scratch)
            os.rmdir(model_dir)
        except FileNotFoundError:
            pass
        os.makedirs(new_scratch)
        if new_scratch != model_dir:
            os.mkdir(model_dir)

        with open(os.path.join(constants.OUTPUT_DIR, 'TensorFlow %s Training Accuracy.txt' % category_name), 'w') as f:

            f.write('Phenotype category: %s\n' % category_name)

            # can adjust activation function (sigmoid works okay)
            # very easy to overfit though; not enough data to train more sophisticated NN
            # one layer seems to work just as well as more complex models
            # less prone to latch-up: predicting all false or true for a given label (generall false for rare labels)

            classifier = tf.estimator.DNNClassifier(feature_columns=feature_columns,
                                                    hidden_units=[32],
                                                    activation_fn=tf.sigmoid,
                                                    n_classes=2,
                                                    model_dir=new_scratch)

            classifier.train(input_fn=train_input, steps=500)
            f.write(str(classifier.evaluate(input_fn=test_input, steps=50)) + '\n')


def load_tensorflow_model_paths(model_dir, categories):
    category_to_model = dict()

    for cat in categories:
        category_to_model[cat] = os.path.join(model_dir, cat)

    return category_to_model


def cross_categorize_dataset(cat_to_model_dict, feature_type, category_exclusions):
    """

    Attempts to evaluate correlations between NN predicts statistically.

    :param cat_to_model_dict: 
    :param feature_type: 
    :param category_exclusions: 
    :return: 
    """

    import itertools

    def heatmap_constructor(x_y_heatmap):

        heatmap = defaultdict(dict)

        for x in label_names.values():
            for y in label_names.values():
                heatmap[x][y] = 0

        for j in x_y_heatmap:

            for (p1, p2) in itertools.product(x_y_heatmap[j], x_y_heatmap[j]):
                if p1 == p2:
                    continue
                else:
                    heatmap[p1][p2] += 1

        return heatmap

    def matrix_constructor(heatmap, id_pairs, norm_factor):

        matrix = []

        for (phenotype_x, _) in id_pairs:
            temp = []
            for (phenotype_y, _) in id_pairs:
                # normalize counts by total genotypes
                temp.append(float(heatmap[phenotype_x][phenotype_y]) / float(norm_factor))
            matrix.append(temp)

        return matrix

    feature_set, feature_dataset, label_data, label_names = generate_tf_features_labels(feature_type,
                                                                                        category_exclusions)

    feature_columns = [tf.contrib.layers.real_valued_column('test',
                                                            dimension=len(feature_set),
                                                            dtype=tf.int32)]

    predicted_phenotypes = defaultdict(list)
    actual_phenotypes = defaultdict(list)
    correlation_record_dict = defaultdict(list)

    def input_fn():

        features = feature_dataset
        batch_size = 1

        dataset = {'test': tf.convert_to_tensor(features, dtype=tf.int32)}
        dataset = tf.data.Dataset.from_tensor_slices(dataset)
        dataset = dataset.batch(batch_size)

        return dataset

    for i in label_names:

        category = label_names[i]

        if category not in cat_to_model_dict:
            continue

        model_path = cat_to_model_dict[category]

        with tf.Session() as session:

            # NN hidden units must match that used for training
            classifier = tf.estimator.DNNClassifier(feature_columns=feature_columns,
                                                    hidden_units=[32],
                                                    model_dir=model_path)

            result = classifier.predict(input_fn=input_fn)

            counter = 0

            for x in result:

                if label_data[i][counter] == 1:
                    actual_phenotypes[counter].append(category)
                if int(x['classes'][0]) == 1:
                    predicted_phenotypes[counter].append(category)
                    correlation_record_dict[category].append(1)
                else:
                    correlation_record_dict[category].append(0)
                counter += 1

    predicted_heatmap = heatmap_constructor(predicted_phenotypes)
    actual_heatmap = heatmap_constructor(actual_phenotypes)

    rows = label_names.values()
    converted_identifiers = sql_interface.convert_full_identifiers_to_abbreviations('phenotype',
                                                                                    rows)
    id_pairs = sorted([(x, y) for x, y in zip(rows, converted_identifiers)],
                      key=lambda k: k[1],
                      reverse=True)

    fish_exact_dict = defaultdict(dict)

    # test enrichment of phenotype pairs versus actual labels in resistome
    for (px, _) in id_pairs:
        for (py, _) in id_pairs:

            if px == py:
                fish_exact_dict[px][py] = 1.0

            sum_px = sum([predicted_heatmap[px][y] for y in predicted_heatmap[px]])
            sum_py = sum([predicted_heatmap[py][y] for y in predicted_heatmap[py]])
            sum_ax = sum([actual_heatmap[px][y] for y in actual_heatmap[px]])
            sum_ay = sum([actual_heatmap[py][y] for y in actual_heatmap[py]])

            # cols: predicted / actual
            # rows: same strain / not same strain
            contingency_table = [[predicted_heatmap[px][py], actual_heatmap[px][py]],
                                 [sum_px + sum_py - predicted_heatmap[px][py],
                                  sum_ax + sum_ay - actual_heatmap[px][py]]]

            (oddsratio, p) = scipy.stats.fisher_exact(contingency_table)

            if numpy.isnan(p):
                fish_exact_dict[px][py] = 0
            elif p == 0:
                fish_exact_dict[px][py] = 25
            else:
                fish_exact_dict[px][py] = min(abs(numpy.log10(p)), 25)

    predicted_matrix = matrix_constructor(fish_exact_dict, id_pairs, 1.0)
    # visualize prediction correlation frequencies (not statistical)
    glib.generate_heatmap(predicted_matrix,
                          [k[1] for k in id_pairs],
                          'Phenotype Associations', os.path.join(constants.OUTPUT_DIR,
                                                                 'Tensorflow NN Phenotype Correlation Counts.pdf'),
                          cmap='Greens')

    predicted_matrix = []

    for (phenotype_x, _) in id_pairs:
        temp = []
        for (phenotype_y, _) in id_pairs:

            if phenotype_x == phenotype_y:
                R = 0.0
            else:
                print (phenotype_x, phenotype_y)
                (R, p) = scipy.stats.spearmanr(correlation_record_dict[phenotype_x],
                                               correlation_record_dict[phenotype_y])
            temp.append(R ** 2)
        predicted_matrix.append(temp)

    glib.generate_heatmap(predicted_matrix,
                          [k[1] for k in id_pairs],
                          'Phenotype Associations', os.path.join(constants.OUTPUT_DIR,
                                                                 'Tensorflow NN Phenotype Correlation Spearman.pdf'),
                          cmap='viridis')


if __name__ == '__main__':

    retrain = False

    if retrain:
        train_tensorflow_model(feature_type='go',
                               scratch_dir=os.path.join(constants.INPUT_DIR,
                                                        'neural_network_models'))

    pclasses = sql_interface.get_phenotype_classes()
    converted_pclasses = sql_interface.convert_full_identifiers_to_abbreviations(abbrev_type='phenotype',
                                                                                 keys=pclasses)

    output_data = []

    for pclass, conv_pclass in zip(pclasses, converted_pclasses):

        fname = 'TensorFlow %s Training Accuracy.txt' % pclass

        with open(os.path.join(constants.OUTPUT_DIR, fname)) as f:
            f.readline()
            output_accuracy_data = f.readline().strip()
            output_accuracy_data = output_accuracy_data.replace('{', '').replace('}', '')

            try:
                accuracy = output_accuracy_data.split(',')[0].split(':')[1]
                output_data.append((conv_pclass, float(accuracy)))
            except:
                pass

    output_data = sorted(output_data, key=lambda x: x[0])

    glib.bargraph([x[0] for x in output_data],
                  [x[1] for x in output_data],
                  'Phenotype Class',
                  'Prediction Accuracy',
                  os.path.join(constants.OUTPUT_DIR, 'TensorFlow Prediction Accuracy.pdf'),
                  rotation='vertical')

    # run_analysis()

    cat_model_dict = load_tensorflow_model_paths(os.path.join(constants.INPUT_DIR,
                                                              'neural_network_models'),
                                                 categories={'solvents_biofuels',
                                                             'radiation',
                                                             'oxidative',
                                                             'other',
                                                             'osmotic_stress',
                                                             'organic_acid (neutral)',
                                                             'nutrient_limitation',
                                                             'mutagens',
                                                             'metal_ions',
                                                             'low_temperature',
                                                             'low_ph',
                                                             'high_temperature',
                                                             'high_ph',
                                                             'high_gravity',
                                                             'general_growth',
                                                             'general_antimetabolite',
                                                             'furans',
                                                             'detergents',
                                                             'aa_antimetabolite',
                                                             'antichemo'})

    cross_categorize_dataset(cat_model_dict,
                             feature_type='go',
                             category_exclusions=set(['quaternary_ammonia']))
