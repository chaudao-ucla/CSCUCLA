/* 
 * ALCTHelperFunctions.cpp
 * 
 * Created on: 7 August 2019
 *      Author: Chau Dao 
 * 
 * Adapted from: https://github.com/cms-sw/cmssw/blob/master/L1Trigger/CSCTriggerPrimitives/src/CSCAnodeLCTProcessor.cc
 */


#include "../include/ALCTHelperFunctions.h"
#include "../include/CSCClasses.h"
#include "../include/CSCConstants.h"

std::vector<int> pulse_to_vec(const unsigned int pulse)
{
    std::vector<int> tbins;
    uint32_t tbit = pulse;
    uint32_t one = 1;
    for (int i = 0; i < 32; i++) 
    {
        if (tbit & one) tbins.push_back(i);
        tbit = tbit >> 1;
        if (tbit == 0) break;
    }
  return tbins;
}

void print_pulse(unsigned int pulse)
{   
    uint32_t tbit = pulse;
    uint32_t one = 1;
    std::cout << "Begin:";
    int prev_state = 0;
    for (int i = 0; i < 32; i++)
    {
        if (tbit & one)
        {   
            if (prev_state == 0) {std::cout << "|-"; prev_state = 1;}
            else if (prev_state == 1) std::cout << "-";
        }
        else
        {
            if (prev_state == 0) std::cout << "_"; 
            else if (prev_state == 1) {std::cout << "|_"; prev_state = 0;}
        }
        tbit = tbit >> 1; 
    }
    std::vector<int> vec_pulse = pulse_to_vec(pulse);
    for (int i=0; i<vec_pulse.size(); i++)
    {
        std::cout << vec_pulse.at(i) << " ";
    }
    std::cout << std::endl;
    std::cout << std::endl;
}

unsigned int extend_time(const unsigned int pulse, const int p_ext)
{
    uint32_t tbit = pulse;
    uint32_t one = 1; 
    for (int i = 0; i < 16; i++) 
    {
        if (tbit & one)
        {
            for (int j = 0; j<p_ext-1; j++)
            {
                one <<= 1; 
                tbit|=one; 
                i++;
                if (i>=16) return tbit;
            }
        }
        one <<= 1;
    }
  return tbit;
}

bool preTrigger(std::vector<ALCT_ChamberHits*> &chamber_list, 
                ALCTConfig &config,
                ALCTCandidate &cand)
{
    int kwg = cand.get_kwg();
    int start_bx = config.get_start_bx();
    int i_pattern = cand.get_pattern();
    int pattern_mask[N_ALCT_PATTERNS][MAX_WIRES_IN_PATTERN];

    for (int i_patt = 0; i_patt < N_ALCT_PATTERNS; i_patt++) 
    {
        for (int i_wire = 0; i_wire < MAX_WIRES_IN_PATTERN; i_wire++) 
        {
            pattern_mask[i_patt][i_wire] = pattern_mask_open[i_patt][i_wire];
            if (((chamber_list.at(0))->_ring== 1 || (chamber_list.at(0))->_ring == 4))
                if (config.narrow_mask_flag())
                    pattern_mask[i_patt][i_wire] = pattern_mask_r1[i_patt][i_wire];
                else 
                    pattern_mask[i_patt][i_wire] = pattern_mask_open[i_patt][i_wire];
        }   
    }

    int layers_hit;
    bool hit_layer[NLAYERS];

    const int nplanes_hit_pretrig_acc = 
        (config.get_nplanes_accel_pretrig() != 0) ? 
            config.get_nplanes_accel_pretrig() : 
            config.get_nplanes_hit_pretrig();
    const int pretrig_thresh[N_ALCT_PATTERNS] = 
    {
        nplanes_hit_pretrig_acc, 
        config.get_nplanes_hit_pretrig(), 
        config.get_nplanes_hit_pretrig()
    };

    unsigned int stop_bx = config.get_fifo_tbins() - config.get_drift_delay();
    for (int i = start_bx; i<stop_bx; i++)
    {
        ALCT_ChamberHits* chamber = chamber_list.at(i);
        if (chamber->get_nhits()<pretrig_thresh[i_pattern]) continue; 
        int MESelect = (chamber->_station <= 2) ? 0 : 1;
        for (int i_lay = 0; i_lay<NLAYERS; i_lay++)
            hit_layer[i_lay] = false;
        layers_hit = 0;
        for (int i_wire = 0; i_wire < MAX_WIRES_IN_PATTERN; i_wire++)
        {
            if (!pattern_mask[i_pattern][i_wire]) continue; 
            int this_layer = pattern_envelope[0][i_wire];
            int this_wire = pattern_envelope[1 + MESelect][i_wire] + kwg;
            if (this_wire < 0 || this_wire >= chamber->get_maxWi()) continue;
            if (chamber->_hits[this_wire][this_layer] && !hit_layer[this_layer])
            {
                hit_layer[this_layer] = true;
                layers_hit++;
            }
            if (layers_hit >= pretrig_thresh[i_pattern])
            {
                cand.set_first_bx(i);
                return true;
            }
        }
    }
    cand.nix();
    return false;
}

void preTrigger(std::vector<ALCT_ChamberHits*> &chamber_list, 
                ALCTConfig &config,
                ALCTCandidate * &head)
{   
    if (head==NULL) return;
    ALCTCandidate * cand = head; 
    while (cand!= NULL)
    {
        ALCTCandidate * temp = cand->next; 
        bool trigger = preTrigger(chamber_list, config, *cand);
        if (!trigger && cand == head) head = temp;
        cand = temp;
    }
    return; 
}

bool patternDetection(  const std::vector<ALCT_ChamberHits*> &chamber_list, 
                        const ALCTConfig &config,
                        ALCTCandidate &cand)
{
    int key_wire = cand.get_kwg();
    int i_pattern = cand.get_pattern(); 
    int pattern_mask[N_ALCT_PATTERNS][MAX_WIRES_IN_PATTERN];

    for (int i_patt = 0; i_patt < N_ALCT_PATTERNS; i_patt++) 
    {
        for (int i_wire = 0; i_wire < MAX_WIRES_IN_PATTERN; i_wire++) 
        {
            pattern_mask[i_patt][i_wire] = pattern_mask_open[i_patt][i_wire];
            if (((chamber_list.at(0))->_ring== 1 || (chamber_list.at(0))->_ring == 4))
                if (config.narrow_mask_flag())
                    pattern_mask[i_patt][i_wire] = pattern_mask_r1[i_patt][i_wire];
                else 
                    pattern_mask[i_patt][i_wire] = pattern_mask_open[i_patt][i_wire];
        }   
    }

    bool trigger = false;
    bool hit_layer[NLAYERS];
    unsigned int temp_quality;

    int this_layer, this_wire, delta_wire;

    const int nplanes_hit_pattern_acc = 
        (config.get_nplanes_accel_pattern() != 0) ? 
            config.get_nplanes_accel_pattern() : 
            config.get_nplanes_hit_pattern();

    const int pattern_thresh[N_ALCT_PATTERNS] = 
    {
        nplanes_hit_pattern_acc, 
        config.get_nplanes_hit_pattern(), 
        config.get_nplanes_hit_pattern()
    };

    ALCT_ChamberHits* chamber = chamber_list.at(0);
    int MESelect = (chamber->_station <= 2) ? 0 : 1;
    temp_quality = 0;
    for (int i_layer = 0; i_layer < NLAYERS; i_layer++)
        hit_layer[i_layer] = false; 
    
    double times_sum = 0.;
    double num_pattern_hits = 0.;
    std::multiset<int> mset_for_median;
    mset_for_median.clear();

    for (int i_wire = 0; i_wire < MAX_WIRES_IN_PATTERN; i_wire++)
    {
        if (pattern_mask[i_pattern][i_wire])
        {
            this_layer = pattern_envelope[0][i_wire];
            delta_wire = pattern_envelope[1 + MESelect][i_wire];
            this_wire = delta_wire + key_wire;

            if (this_wire<0 || this_wire>= chamber->get_maxWi()) continue;
            chamber = chamber_list.at(cand.get_first_bx()+config.get_drift_delay());
            if (chamber->_hits[this_wire][this_layer])
            {
                if (!hit_layer[this_layer])
                {
                    hit_layer[this_layer] = true;
                    temp_quality++;
                }
                if (abs(delta_wire)<2)
                {
                    ALCT_ChamberHits* temp = chamber; 
                    int first_bx_layer = cand.get_first_bx() + config.get_drift_delay();
                    for (unsigned int dbx = 0; dbx < config.get_hit_persist(); dbx++)
                    {
                        if (chamber->prev == NULL) break;
                        temp = temp->prev;
                        if (temp->_hits[this_wire][this_layer]) first_bx_layer--;
                        else break;    
                    }
                    times_sum += (double) first_bx_layer;
                    num_pattern_hits += 1.; 
                    mset_for_median.insert(first_bx_layer);
                }
            }
        }
    }
    const int sz = mset_for_median.size();
    if (sz > 0) 
    {
        std::multiset<int>::iterator im = mset_for_median.begin();
        if (sz > 1) std::advance(im, sz / 2 - 1);
        if (sz == 1) cand.set_first_bx_corr(*im);
        else if ((sz % 2) == 1) cand.set_first_bx_corr(*(++im));
        else cand.set_first_bx_corr(((*im) + (*(++im))) / 2);
    }
    if (temp_quality >= pattern_thresh[i_pattern]) 
    {
        trigger = true;
        temp_quality = getTempALCTQuality(temp_quality);

        if (i_pattern == 0) 
        {
            // Accelerator pattern
            cand.set_quality(temp_quality);
            cand.set_pattern(i_pattern);
        }
        else 
        {
            // Only one collision pattern (of the best quality) is reported
            if (static_cast<int>(temp_quality) > cand.get_quality()) 
            {
                cand.set_quality(temp_quality); //real quality
                cand.set_pattern(i_pattern); // pattern, left or right
            }
        }
    }
    else cand.nix();
    return trigger;
}

void patternDetection(  std::vector<ALCT_ChamberHits*> &chamber_list, 
                        ALCTConfig &config,
                        ALCTCandidate * &head)
{
    if (head==NULL) return;
    ALCTCandidate * cand = head; 
    while (cand!= NULL)
    {
        ALCTCandidate * temp = cand->next; 
        bool trigger = patternDetection(chamber_list, config, *cand);
        if (!trigger && cand == head) head = temp;
        cand = temp;
    }
    return; 
}

void ghostBuster(ALCTCandidate* curr)
{
    if (curr == NULL) return;
    ALCTCandidate* temp = curr->next; 
    if (temp == NULL) return; 
    if (curr->get_quality())
    {
        if (temp->get_kwg()-1 == curr->get_kwg())
        {
            int dt = curr->get_first_bx() - temp->get_first_bx();
            if (dt == 0)
            {
                if (curr->get_quality()>=temp->get_quality()) temp->flag();
                else curr->flag(); 
            }
            else if (dt<=ghost_cancel && curr->get_quality() > temp->get_quality())
                temp->flag();
            else if (dt<=ghost_cancel && curr->get_quality() < temp->get_quality())
                curr->flag();
        }
    }
    else curr->flag(); 
    ghostBuster(temp); 
}

void clean(ALCTCandidate* &curr)
{
    if (curr == NULL) return; 
    if (curr->prev == NULL) // if this is the current head of list
    {
        if (!curr->isValid()) // if the head of the list is invalid
        {
            ALCTCandidate * temp = curr; 
            curr = curr->next;
            temp->nix();
            clean(curr);
        }
        else
        {
            clean(curr->next);
        }
    }
    else // current element is not the head and is not NULL
    {
        ALCTCandidate* temp = curr->next;
        if (!curr->isValid()) curr->nix();
        clean(temp);
    }
}

void head_to_vec(ALCTCandidate* curr, std::vector<ALCTCandidate*> &cand_vec)
{
    if (curr == NULL) return; 
    cand_vec.push_back(curr);
    head_to_vec(curr->next, cand_vec);
}

int getTempALCTQuality(int temp_quality)
{
    return (temp_quality > 3) ? temp_quality - 3 : 0; 
}

void wipe(std::vector<ALCTCandidate*> candvec)
{
    for (int i=0; i<candvec.size(); i++)
    {
        delete candvec.at(i);
    }
}

void wipe(std::vector<ALCT_ChamberHits*>cvec)
{
    for (int i=0; i<cvec.size(); i++)
    {
        delete cvec.at(i);
    }
}

void bestTrackSelector()
{

}



