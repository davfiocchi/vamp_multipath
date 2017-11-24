// This is a skeleton file for use in creating your own plugin
// libraries.  Replace MyPlugin and myPlugin throughout with the name
//   of your first plugin class, and fill in the gaps as appropriate.


#include "MyPlugin.h"

MyPlugin::MyPlugin(float inputSampleRate) :
	Plugin(inputSampleRate),
	m_blockSize(0),
	m_sampleRate(inputSampleRate),
	m_stepSize(0),
	durationInputFile(0),
	n(0)
{
}

MyPlugin::~MyPlugin()
{
}

string
MyPlugin::getIdentifier() const
{
	return "beattrackerMP";
}

string
MyPlugin::getName() const
{
	return "Multipath Beat Tracker";
}

string
MyPlugin::getDescription() const
{
	// Return something helpful here!
	return "Beat Tracker PlugIn";
}

string
MyPlugin::getMaker() const
{
	// Your name here
	return "Davide, Mattia, Alessandro";
}

int
MyPlugin::getPluginVersion() const
{
	// Increment this each time you release a version that behaves
	// differently from the previous one
	return 1.0;
}

string
MyPlugin::getCopyright() const
{
	// This function is not ideally named.  It does not necessarily
	// need to say who made the plugin -- getMaker does that -- but it
	// should indicate the terms under which it is distributed.  For
	// example, "Copyright (year). All Rights Reserved", or "GPL"
	return "Polimi";
}

MyPlugin::InputDomain
MyPlugin::getInputDomain() const
{
	return FrequencyDomain;
}

size_t
MyPlugin::getPreferredBlockSize() const
{
	return 1024; // 0 means "I can handle any block size"
}

size_t
MyPlugin::getPreferredStepSize() const
{
	return 512; // 0 means "anything sensible"; in practice this
				// means the same as the block size for TimeDomain
				// plugins, or half of it for FrequencyDomain plugins
}

size_t
MyPlugin::getMinChannelCount() const
{
	return 1;
}

size_t
MyPlugin::getMaxChannelCount() const
{
	return 2;
}

MyPlugin::ParameterList
MyPlugin::getParameterDescriptors() const
{
	ParameterList list;

	// If the plugin has no adjustable parameters, return an empty
	// list here (and there's no need to provide implementations of
	// getParameter and setParameter in that case either).

	// Note that it is your responsibility to make sure the parameters
	// start off having their default values (e.g. in the constructor
	// above).  The host needs to know the default value so it can do
	// things like provide a "reset to default" function, but it will
	// not explicitly set your parameters to their defaults for you if
	// they have not changed in the mean time.

	ParameterDescriptor ntrackers;
	ntrackers.identifier = "ntrackers";
	ntrackers.name = "Number of trackers";
	ntrackers.description = "The number of trackers that the algorithm should use";
	ntrackers.unit = "";
	ntrackers.minValue = 2;
	ntrackers.maxValue = 40;
	ntrackers.defaultValue = 18;
	ntrackers.isQuantized = true;
	ntrackers.quantizeStep = 1.0f; //just integer values
	list.push_back(ntrackers);

	ParameterDescriptor stability;
	stability.identifier = "stability";
	stability.name = "Trackers stability";
	stability.description = "The stability of the trackers";
	stability.unit = "";
	stability.minValue = 10;
	stability.maxValue = 100;
	stability.defaultValue = 40;
	stability.isQuantized = true;
	stability.quantizeStep = 0.5f;
	list.push_back(stability);

	return list;
}

float
MyPlugin::getParameter(string identifier) const
{
	if (identifier == "ntrackers") return m_ntrackers;
	else if (identifier == "stability") return m_stability;
	return 10;
}

void
MyPlugin::setParameter(string identifier, float value)
{
	if (identifier == "ntrackers") {
		if ((value >= 1) && (value <= 40))
			m_ntrackers = value;
	}
	else if (identifier == "stability") {
		if ((value >= 10) && (value <= 100))
			m_stability = value;
	}
}

MyPlugin::ProgramList
MyPlugin::getPrograms() const
{
	ProgramList list;

	// If you have no programs, return an empty list (or simply don't
	// implement this function or getCurrentProgram/selectProgram)

	return list;
}

string
MyPlugin::getCurrentProgram() const
{
	return ""; // no programs
}

void
MyPlugin::selectProgram(string name)
{
}

MyPlugin::OutputList
MyPlugin::getOutputDescriptors() const
{
	OutputList list;

	// See OutputDescriptor documentation for the possibilities here.
	// Every plugin must have at least one output.

	OutputDescriptor beats;
	beats.identifier = "beats";
	beats.name = "Beats";
	beats.description = "The beat of the input file";
	beats.unit = "";
	beats.hasFixedBinCount = true;
	beats.binCount = 0;
	beats.hasKnownExtents = false;
	beats.isQuantized = false;
	beats.sampleType = OutputDescriptor::VariableSampleRate;
	beats.sampleRate = (m_inputSampleRate / m_stepSize);
	beats.hasDuration = false;

	list.push_back(beats);

	return list;
}

bool
MyPlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
	if (channels < getMinChannelCount() ||
		channels > getMaxChannelCount()) return false;

	// Real initialisation work goes here!
	m_blockSize = blockSize;
	m_stepSize = stepSize;
	frame_2 = zeros<cx_rowvec>(m_blockSize);
	frame_1 = zeros<cx_rowvec>(m_blockSize);
	envelopeArray = zeros<rowvec>(0);
	origin = Vamp::RealTime::zeroTime;
	return true;
}

void
MyPlugin::reset()
{
	// Clear buffers, reset stored values, etc
}

MyPlugin::FeatureSet
MyPlugin::process(const float *const *inputBuffers, Vamp::RealTime timestamp)
{
	Feature f;
	FeatureSet fs;
	double env;
	real = zeros<rowvec>(m_blockSize);
	img = zeros<rowvec>(m_blockSize);
	int k = 0;

	for (size_t i = 0; i<m_blockSize + 2; i += 2)
	{
		real(k) = inputBuffers[0][i];
		img(k) = inputBuffers[0][i + 1];
		++k;
	}

	// durationInputFile = durationInputFile + 1;

	frame = cx_rowvec(real, img);
	env = getODFValue(frame, frame_2);
	frame_2 = frame_1;
	frame_1 = frame;

	storeODF.push_back(env);

	return fs;
}

MyPlugin::FeatureSet
MyPlugin::getRemainingFeatures()
{
	FeatureSet fs;
	rowvec ibi = getBeatPeriod(storeODF);
	rowvec odfNorm = normalizeODF(storeODF, ibi.max());

	size_t n = odfNorm.n_elem;
	int m = 3;	//indexes, ODF and IBI
	int* out_path;	//where to store the beat index
	int path_len;	//length of out_path
	int* peak_beat = new int[n];

	//calculate first and last index for the valid range
	int firstIdx = 0;
	int lastIdx;
	for (lastIdx = odfNorm.n_elem - 1; (lastIdx > 0) && (odfNorm(lastIdx)<0); lastIdx--);
	n = lastIdx - firstIdx;

	//build the array to pass to the Multipath
	float *mixedArray = new float[3 * n];
	for (int i = 0; i<m; i++) {
		for (int j = firstIdx; j<lastIdx; j++) {
			if (i == 0) {
				mixedArray[i*n + j] = j;
			}
			else if (i == 1) {
				mixedArray[i*n + j] = odfNorm(j);
			}
			else if (i == 2) {
				mixedArray[i*n + j] = ibi(j);
			}
			peak_beat[j] = 0;
		}
	}

	PathFinderMP *pathFinder = new PathFinderMP(m_ntrackers, m_stability);
	pathFinder->setPtsList(mixedArray, m, n);
	pathFinder->findPath(&out_path, &path_len);

	for (int i = 0; i<path_len; i++) {
		int beat_index = out_path[i];
		peak_beat[beat_index] = 150;
	}
	cout << "PATH LENGTH: " << path_len << endl;
	
	//return beats
	for (size_t i = firstIdx; i<lastIdx; i++)
	{
		if (peak_beat[i] > 0) { //if we have a beat return a feature
			Feature f;
			size_t frame = i*m_stepSize;
			f.hasTimestamp = true;
			f.timestamp = origin + Vamp::RealTime::frame2RealTime(frame, lrint(m_sampleRate));
			fs[0].push_back(f);
		}
	}

	return fs;
}

