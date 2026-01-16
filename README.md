# An orbitofrontal microcircuit for approach-avoidance decision-making

Dataset DOI: [10.5061/dryad.kh18932k7](10.5061/dryad.kh18932k7)

Code reproducing main figures. All electrophysiological and behavioral data is contained in subject_data.mat, available at Dryad (https://doi.org/10.5061/dryad.kh18932k7) following completion of the peer review process.

There are example outputs of the decoding and posterior generating scripts as these take a while to run, see "posteriors_LOO" folder for these examples.

## Description of the data and file structure

### Files and variables

#### File: subject_data.mat

**Description:** 

##### Variables

* subject is the "master" variable. All main figure analyses from the manuscript can be reproduced using the "subject" variable. Below I explain the subfields. Use subject(n).subfield1 to obtain subfield1 information for subject n
* STORAGE OF TRIAL INFORMATION AND BEHAVIOR:
  * decision: the decision made on a particular trial. ordered sequentially from 1-n_trials. 1 = approach. 0 = avoidance.
  * The below variables (all subfields of "subject") are the key to understanding the task the subject played, and overall how they played it:
    *     reward_trial_type: for trial type n, reward_trial_type(n) is the reward offer for that trial type (# lit treasures)

          punishment_trial_type: for trial type n, punishment_trial_type(n) is the punishment offer for that trial type (# lit bombs)

          conflict_trial_type: for trial type n, conflict_trial_type(n) is the behavioral conflict the subject exhibited

          trial_type_trial_type: the numbered trial type.

          p_approach_trial_type: for trial type n, p_approach_trial_type(n) is the subject's overall probability of approach
  * The below variables give the trial-by-trial log of the subject's behavior. Each will be 1-ntrials long.
    *         p_approach_trial: subject's overall approach probability for the particular trial type presented in the nth trial 

          conflict_trial: subject's overall behavioral conflict for the particular trial type presented in the nth trial

          reward_trial: reward offer on a particular trial

          punish_trial: punishment offer on a particular trial

          trial_type_trial: trial type from beginning to end of session
  * trial_type_trial: index of distinct trial type spanning from 1-ntrials. each trial type corresponded to a distinct reward/punishment combo.
* STORAGE OF NEURAL DATA:
  * electrode: information about each electrode for each subject. So, subject(1).electrode(1) provides coordinates (lat_coor/med_coor/olf_coor/trans_coor) relative to lateral orbital sulcus, medial orbital sulcus, olfactory sulcus, and transverse orbital sulcus for electrode 1 in subject 1. I suggest adhering to these coordinates rather than "anatomical info" as the coordinates are most accurate and obtained using manual surface labeling by neuroanatomists (see info in Methods for the anatomical method) whereas the anatomical info was not based on these detailed reconstructions.
    * The subfield within electrode is trigger: 1 = align to trial onset. 2 = align to decision (button press).
      * The subfield within trigger is high_gamma_mat, which is the neural data (high frequency activity). So, for instance, subject(1).electrode(1).trigger(1).high_gamma_mat gives a high_gamma_mat: [220×12000 double] where each row is a trial from 1-ntrials and each column is a timepoint (in milliseconds) from that trial, aligned to the trigger of interest. So, for example, this would be high frequency activity for the electrode 1 in subject 1 spanning from 6000ms prior to trial onset to 6000ms after trial onset.




Methods for data collection and preprocessing (this text is reproduced from the submitted manuscript that is currently under revision):

Participants

Data were obtained from 6 patients (Table S1) who had stereotactically implanted intracranial depth electrodes (Ad-Tech or Dixi, 3-5mm inter-electrode spacing) placed at the University of California, San Francisco (UCSF, n = 5 subjects) or at Washington University in St. Louis (WUSTL, n = 1 subject). Every patient with stereo-electroencephalography implants in OFC at UCSF between June 2023-September 2024, and at WUSTL between August 2024-October 2024, was included in the dataset. Each subject was treated as a biological replicate. No subject randomization was performed as all subjects performed the same task. Electrode placement was exclusively guided by clinical needs. Approval for the study was granted by the institutional review boards of the UCSF, UC Berkeley and WUSTL. Written informed consent was obtained by all subjects prior to testing.

 

Behavior and task

Participants each completed one run of the approach-avoidance task (playable at neurogame.ucsf.edu). The task was coded in the cross-platform game engine Unity (Unity Technologies, San Francisco, CA) by an independent gaming studio (Cat-astrophe Games, Wrocław, Poland).

 

The task took approximately 20 minutes for subjects to complete (192-220 trials). Participants had 6 seconds to indicate approach or avoidance decisions with a button press, or else a single ruby would be subtracted and the trial would automatically end. Trial types were randomly interleaved and pseudorandomized across blocks of 40 trials to avoid clumping over a particular trial type during the experiment. The intertrial interval was drawn from a Poisson distribution (median = 985ms, S.D = 249.41ms). For Subjects 1, 2, 3, 4, and 6 (experiment run at UCSF) Black Box ToolKit mBBTK v2/event marking model Elite (Black Box ToolKit, Sheffield, United Kingdom) was used to record task-related events. Trial starts were recorded using a photodiode. Button presses were recorded with the BlackBox ToolKit USB response pad. The precise timing of outcomes (bomb explosions, treasure chests opening) was accompanied by auditory output recorded with a microphone. These were conveyed as analog inputs to the recording system and used to synchronize behavioral and electrophysiological data. For Subject 5 (experiment run at WUSTL), the Brain Computer Interface 2000 (BCI2000) was used to record button presses from the laptop keyboard pad and trial starts from a photodiode. The participants viewed the task on a screen placed in front of them and used the BBTK USB response pad (or laptop, in the case of Subject 5), to indicate approach versus avoid. Subjects 1, 4, 5, and 6 played a 220-trial version of the task. This included reward magnitudes of 0, 2, 4, and 7; and punishment magnitudes of 0, 2, and 4. Every possible permutation of these reward and punishment magnitudes constituted a trial type. Subjects 2 and 3 played the 192-trial version of the task. This included the same range of reward and punishment offers with more granularity: reward magnitudes of 0,1,2,3,4,5,6,7 and punishment magnitudes of 0,1,2,3,4 were presented.



Manual definition of orbitofrontal sulci

Each structural T1-weighted magnetic resonance imaging (MRI) scan was first processed through FreeSurfer v6.0.0 to create 3D pial and inflated cortical reconstructions (*31*, *32*). Next, all sulci were defined on these cortical reconstructions using the FreeSurfer curvature metric to differentiate between sulcal and gyral components (*31*–*33*). Trained neuroanatomists (E.H.W. and K.S.W.) manually identified OFC sulci in each hemisphere (N = 6; 12 hemispheres) using FreeSurfer's tksurfer tools and guided by the most up-to-date definitions of OFC sulcal organization (*34*, *35*). The OFC sulci of interest in the present study were as follows: (1) olfactory sulcus (OlfS), (2) transverse olfactory sulcus (TOlfS), (3) transverse orbital sulcus (TOS), (4) anterior section of the medial orbital sulcus (MOS-a), (5) posterior section of the medial orbital sulcus (MOS-p), (6) anterior section of the lateral orbital sulcus (LOS-a), (7) posterior section of the lateral orbital sulcus (LOS-p), (8) intermediate orbital sulcus (IOS), (9) posterior orbital sulcus (POS), and (10) sulcus fragmentosus (SF). Prior research has shown that the OlfS, TOlfS, TOS, MOS-a, MOS-p, LOS-a, and LOS-p are present in all individuals, whereas the IOS, POS, and SF are more variably present, consisting of a variable number of components (ranging from 0-4) (*15*, *34*). The seven consistent sulci were present in all patients. Example sulcal labels for constant sulci MOS (includes MOS-a and MOS-p), LOS (includes LOS-a and LOS-p), TOS, and OlfS are shown in Fig. S2A. The incidence rates of the three variable sulci in each patient are presented in Table S2.

 

Electrode localization

Electrodes were localized based on co-registration of 1mm slice thickness T1 weighted MRI and 1mm slice thickness bone window computed tomography (CT) scans. Following co-registration in Freesurfer, electrode locations were marked using the region of interest (ROI) tool. The centroid coordinates, as well as coordinates from manually defined OFC sulcal labels, were exported in patient space for further analysis in MATLAB (Mathworks, Natick, Massachusetts). We created a standardized coordinate system across subjects based on the constant sulci MOS (combining MOS-a and MOS-p) and TOS(*16*) (Fig. S2B). The anterior/posterior coordinate was measured in relation to the TOS. For each electrode, we calculated the distance between the anterior/posterior electrode coordinate and the anterior/posterior coordinate of the centroid of the points contained in the TOS label that lay within the same medial/lateral dimension as the electrode. If the electrode lay beyond the medial/lateral dimension of the TOS, we measured the anterior/posterior distance from the closest TOS medial/lateral slice containing the TOS. The medial/lateral coordinate was measured in relation to the MOS. For each electrode, we calculated the distance between the medial/lateral electrode coordinate and the medial/lateral coordinate of the centroid of the points contained in the MOS label that lay within the same anterior/posterior dimension as the electrode. For anatomic localization, electrodes were excluded if they lay more than 1mm medial to the centroid of the olfactory sulcus (placing the electrode diameter fully within the gyrus rectus); electrodes were considered to be in the medial OFC if they were otherwise medial to MOS; in the anterior OFC if they were lateral to the MOS, anterior to the TOS, and medial to the LOS; in the posterior OFC if they were lateral to the MOS, posterior to the TOS, and medial to the LOS; and in the lateral OFC if they were lateral to the LOS. Electrodes in the inferior frontal gyrus were excluded. Electrode arrays located anterior to TOS (Fig. 2F) were defined as those with the majority of electrodes located anterior to the TOS.

 

Data collection and preprocessing

We recorded from a total of 131 electrode contacts located in the OFC. Intracranial electroencephalography data were acquired using a multichannel pre-amplifier connected to a digital signal processor and digitized at 3-10kHz using either the Tucker-Davis Technologies or Nihon Kohden recording system. There were no seizures recorded during any epochs. Epilepsy patients’ data was manually examined in subsecond epochs during task recording to confirm the absence of epileptiform activity. Furthermore, the OFC was not identified as an epileptiform locus in any of the patients. Data were visually inspected and flat or noisy channels were removed. Bipolar referencing was applied using the nearest neighboring electrode. Electrodes were only included in the dataset if a nearest neighbor for bipolar referencing was available. Following elimination of noisy channels (15 electrodes) and the exclusion of redundant electrodes introduced by bipolar referencing (14 electrodes), 102 electrodes were included in the final dataset. Laplacian and common average referencing were also tested and yielded similar results to those reported in the manuscript. The fieldtrip package in MATLAB was used for the following pre-processing steps. First, the data was bandpass filtered between 0.5 and 200Hz. Next, a notch filter was applied to remove 60Hz noise and harmonics. To compute the amplitude of high frequency activity (HFA – 70-150Hz), the data was first narrowband filtered in 10Hz bands from 70Hz to 150Hz. The Hilbert transform was applied to the filtered data, and the amplitude envelope was computed for each 10Hz band. Each band’s amplitude was normalized to its average amplitude during a 10s pre-task period, to account for the 1/f exponential decay of signal amplitude across frequency bands. The normalized HFA amplitude for each 10Hz band spanning from 70Hz to 150Hz was averaged to provide the final analyzed HFA signal.


## Code/software

Matlab codes to reproduce all main figures in the paper: https://github.com/cstarkweather/OFC_analysis_codes
## Access information

Other publicly accessible locations of the data:

* none

Data was derived from the following sources:

* all collected for this manuscript as stated in Methods



## Human subjects data

There is no personal health information in the deposited data. All data was obtained in accordance with our institutional IRB, and informed consent was obtained from all participants prior to participation.
