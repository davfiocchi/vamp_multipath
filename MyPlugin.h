
 // This is a skeleton file for use in creating your own plugin
// libraries.  Replace MyPlugin and myPlugin throughout with the name
// of your first plugin class, and fill in the gaps as appropriate.


// Remember to use a different guard symbol in each header!
#ifndef MY_PLUGIN_H
#define MY_PLUGIN_H

#include "vamp-sdk/Plugin.h"
#include </usr/local/include/armadillo>
#include <Vector>

#include "getODFValue.h"
#include "IBI.h"//si può rimuovere
#include "peakPicking.h"//si può rimuovere
#include "Multipath/src/PathFinderMP.h"
#include "Multipath/src/SimpleTracker.h"
#include "Multipath/src/TrackersManager.h"
#include "Multipath/src/PtsList.h"
#include "Multipath/src/PathFinder.h"
#include "Multipath/src/BeatPtsList.h"
#include "Multipath/src/kth_order_statistic.h"
#include "BeatPeriod.hpp"


using std::string;
using namespace arma;
using std::vector;

class MyPlugin : public Vamp::Plugin
{
public:
    MyPlugin(float inputSampleRate);
    virtual ~MyPlugin();

    string getIdentifier() const;
    string getName() const;
    string getDescription() const;
    string getMaker() const;
    int getPluginVersion() const;
    string getCopyright() const;

    InputDomain getInputDomain() const;
    size_t getPreferredBlockSize() const;
    size_t getPreferredStepSize() const;
    size_t getMinChannelCount() const;
    size_t getMaxChannelCount() const;

    ParameterList getParameterDescriptors() const;
    float getParameter(string identifier) const;
    void setParameter(string identifier, float value);

    ProgramList getPrograms() const;
    string getCurrentProgram() const;
    void selectProgram(string name);

    OutputList getOutputDescriptors() const;

    bool initialise(size_t channels, size_t stepSize, size_t blockSize);
    void reset();

    FeatureSet process(const float *const *inputBuffers,
                       Vamp::RealTime timestamp);

    FeatureSet getRemainingFeatures();
    
    size_t durationInputFile;

protected:
	// plugin-specific data and methods go here
	size_t m_blockSize;
	size_t m_stepSize;
	float m_sampleRate;

	//Input parameters
	float m_ntrackers = 18;
	float m_stability = 40;

	//Useful variables for the plugin to store information
	cx_rowvec frame;
	cx_rowvec frame_1;
	cx_rowvec frame_2;
	rowvec real;
	rowvec img;
	rowvec peaks;
	rowvec ibi;
	rowvec envelopeArray;
	vector<double> storeODF;
	// FeatureSet fs;
	int j = 0;
	int n;
	vector<Vamp::RealTime> timestamps;
	Vamp::RealTime origin;
};



#endif
