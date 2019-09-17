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
        std::cout << vec_pulse[i] << " ";
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

bool preTrigger(int trig_time,
                std::vector<ALCT_ChamberHits*> &chamber_list, 
                ALCTConfig &config,
                ALCTCandidate &cand)
{
    int kwg = cand.get_kwg();
    int i_pattern = cand.get_pattern();
    int i = trig_time; 
    ALCT_ChamberHits* chamber = chamber_list[i];

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

    if ((chamber->get_nhits()<pretrig_thresh[i_pattern]))
    {
        cand.flag();
        return false; 
    }
    
    int pattern_mask[N_ALCT_PATTERNS][MAX_WIRES_IN_PATTERN];
    for (int i_patt = 0; i_patt < N_ALCT_PATTERNS; i_patt++) 
    {
        for (int i_wire = 0; i_wire < MAX_WIRES_IN_PATTERN; i_wire++) 
        {
            if ((chamber_list[0]->_ring== 1 || chamber_list[0]->_ring == 4))
            {
                if (config.narrow_mask_flag())
                    pattern_mask[i_patt][i_wire] = pattern_mask_r1[i_patt][i_wire];
                else 
                    pattern_mask[i_patt][i_wire] = pattern_mask_open[i_patt][i_wire];
            }
            else 
            {
                pattern_mask[i_patt][i_wire] = pattern_mask_open[i_patt][i_wire];
            }
        }
    }

    int layers_hit;
    bool hit_layer[NLAYERS];

    unsigned int stop_bx = config.get_fifo_tbins() - config.get_drift_delay();
    
    if (!(chamber->get_nhits()<pretrig_thresh[i_pattern]))
    {
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
    cand.flag();
    return false;
}

bool patternDetection(  const std::vector<ALCT_ChamberHits*> &chamber_list, 
                        const ALCTConfig &config,
                        ALCTCandidate &cand)
{
    int key_wire = cand.get_kwg();
    int i_pattern = cand.get_pattern(); 
    int pattern_mask[N_ALCT_PATTERNS][MAX_WIRES_IN_PATTERN];

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

    ALCT_ChamberHits* chamber = chamber_list[cand.get_first_bx()+config.get_drift_delay()];
    if ((chamber->get_nhits()<pattern_thresh[i_pattern]))
    {
        cand.flag();
        return false; 
    }

    for (int i_patt = 0; i_patt < N_ALCT_PATTERNS; i_patt++) 
    {
        for (int i_wire = 0; i_wire < MAX_WIRES_IN_PATTERN; i_wire++) 
        {
            if ((chamber_list[0]->_ring== 1 || chamber_list[0]->_ring == 4))
            {
                if (config.narrow_mask_flag())
                    pattern_mask[i_patt][i_wire] = pattern_mask_r1[i_patt][i_wire];
                else 
                    pattern_mask[i_patt][i_wire] = pattern_mask_open[i_patt][i_wire];
            }
            else 
            {
                pattern_mask[i_patt][i_wire] = pattern_mask_open[i_patt][i_wire];
            }
        }
    }

    bool trigger = false;
    bool hit_layer[NLAYERS];
    unsigned int temp_quality;

    int this_layer, this_wire, delta_wire;

    int MESelect = (chamber->_station <= 2) ? 0 : 1;
    temp_quality = 0;
    for (int i_layer = 0; i_layer < NLAYERS; i_layer++)
        hit_layer[i_layer] = false; 
    
    double times_sum = 0.;
    double num_pattern_hits = 0.;

    for (int i_wire = 0; i_wire < MAX_WIRES_IN_PATTERN; i_wire++)
    {
        if (pattern_mask[i_pattern][i_wire])
        {
            this_layer = pattern_envelope[0][i_wire];
            delta_wire = pattern_envelope[1 + MESelect][i_wire];
            this_wire = delta_wire + key_wire;

            if (this_wire<0 || this_wire>= chamber->get_maxWi()) continue;
            //chamber = chamber_list.at(cand.get_first_bx()+config.get_drift_delay());
            if (chamber->_hits[this_wire][this_layer])
            {
                if (!hit_layer[this_layer])
                {
                    hit_layer[this_layer] = true;
                    temp_quality++;
                }
            }
        }
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
    else cand.flag();
    return trigger;
}

void trig_and_find( std::vector<ALCT_ChamberHits*> &chamber_list, 
                    ALCTConfig &config,
                    std::vector<std::vector<ALCTCandidate*>> &end_vec)
{   
    bool armed = true; 
    for (int j=0; j<chamber_list[0]->get_maxWi(); j++)
    {
        for (int i=0; i<config.get_fifo_tbins()-config.get_drift_delay(); i++)
        {
            ALCTCandidate * curr = end_vec[i][j];
            bool ptrigger = preTrigger(i,chamber_list,config,*curr);
            if (armed && !ptrigger) continue; 
            else if (armed && ptrigger)
            {
                armed = false;
                for (int k=1; k<config.get_drift_delay(); k++)
                {
                    if (i+k<end_vec.size()) end_vec[i+k][j]->flag();
                }
                i+=(config.get_drift_delay());
                bool detect = patternDetection(chamber_list, config, *curr);
                if (!detect) i--;
                //else cout << "detect temp: " << curr << endl << endl;  
            }
            else if (!armed && !ptrigger)
            { 
                armed = true; 
            }
        }
    }
}

void ghostBuster(std::vector<std::vector<ALCTCandidate*>> &end_vec, ALCTConfig &config)
{
    for (int i = 0; i<end_vec.size(); i++)
    {
        std::vector<ALCTCandidate*> wire_vec = end_vec[i];
        for (int j = 0; j<wire_vec.size(); j++)
        {
            ALCTCandidate* this_wire = wire_vec[j];
            if (!this_wire->get_quality()) 
            {
                this_wire->flag();
                continue; 
            }
            if (j!= wire_vec.size()-1)
            {
                ALCTCandidate* next_wire = wire_vec[j+1];
                if (next_wire->get_quality()>this_wire->get_quality())
                {
                    this_wire->flag();
                }
                else next_wire->flag();
            }
            for (int k = 1; k <= config.get_ghost_cancel(); k++)
            {
                if (i-k<0) continue;
                std::vector<ALCTCandidate*> last_time_vec = end_vec[i-k];
                if (j-1 >= 0)
                {
                    ALCTCandidate* last_time = last_time_vec[j-1];
                    bool was_better = last_time->get_quality() > this_wire->get_quality(); 
                    if (last_time->get_quality() && (!config.ghost_flag() || was_better))
                    {
                        this_wire->flag();
                    }
                }
                if (j+1 < wire_vec.size())
                {
                    ALCTCandidate* last_time = last_time_vec[j+1];
                    bool was_better = last_time->get_quality() > this_wire->get_quality(); 
                    if (last_time->get_quality() && (!config.ghost_flag() || was_better))
                    {
                        this_wire->flag();
                    }
                }
            }
        }
    }
}

auto sortRule = [] (ALCTCandidate * cand1, ALCTCandidate * cand2) -> bool
{
    if (cand1->get_quality() == cand2->get_quality()) return (cand1->get_kwg()>cand2->get_kwg());
    else return (cand1->get_quality()>cand2->get_quality()); 
};

auto sortRule1 = [] (ALCTCandidate * cand1, ALCTCandidate * cand2) -> bool
{
    if (cand1->get_quality() == cand2->get_quality()) return (cand1->get_kwg()>cand2->get_kwg());
    else return (cand1->get_quality()>cand2->get_quality()); 
};

auto sortRule2 = [] (ALCTCandidate * cand1, ALCTCandidate * cand2) -> bool
{
    if (cand1->get_quality() == cand2->get_quality()) return (cand1->get_kwg()<cand2->get_kwg());
    else return (cand1->get_quality()>cand2->get_quality()); 
};


void extract(std::vector<std::vector<ALCTCandidate*>> &end_vec, std::vector<ALCTCandidate*> &out_vec)
{
    for (int i = 0; i<end_vec.size(); i++)
    {
        std::vector<ALCTCandidate*> wire_vec = end_vec[i];
        for (int j = 0; j<wire_vec.size(); j++)
        {
            ALCTCandidate * curr = end_vec[i][j];
            if (curr->isValid()) 
            {
                out_vec.push_back(curr);
            }
            else 
            {
                curr->nix();
            }
        }
    }
}

void extract_sort(std::vector<std::vector<ALCTCandidate*>> &end_vec, std::vector<ALCTCandidate*> &out_vec)
{
    for (int i = 0; i<end_vec.size(); i++)
    {
        std::vector<ALCTCandidate*> wire_vec = end_vec[i];
        std::vector<ALCTCandidate*> temp_vec; 
        for (int j = 0; j<wire_vec.size(); j++)
        {
            ALCTCandidate * curr = end_vec[i][j];
            if (curr->isValid()) 
            {
                temp_vec.push_back(curr);
            }
            else 
            {
                curr->nix();
            }
        }
        if (temp_vec.size())
        {
            std::sort(temp_vec.begin(),temp_vec.end(),sortRule);
            for (int k = 0; k<temp_vec.size(); k++)
            {
                if (k==0 || k== 1) out_vec.push_back(temp_vec[k]);
            }
        }
    }
}

void extract_sort_cut(std::vector<std::vector<ALCTCandidate*>> &end_vec, std::vector<ALCTCandidate*> &out_vec)
{
    for (int i = 0; i<end_vec.size(); i++)
    {
        std::vector<ALCTCandidate*> wire_vec = end_vec[i];
        std::vector<ALCTCandidate*> temp_vec; 
        for (int j = 0; j<wire_vec.size(); j++)
        {
            ALCTCandidate * curr = end_vec[i][j];
            if (curr->isValid() && curr->get_first_bx()>=5 && curr->get_first_bx()<=11) 
            {
                temp_vec.push_back(curr);
            }
            else 
            {
                curr->nix();
            }
        }
        if (temp_vec.size())
        {
            int comp_q = -1;
            std::sort(temp_vec.begin(),temp_vec.end(),sortRule1);
            comp_q = temp_vec[0]->get_quality();
            temp_vec[0]->set_tracknumber(1);
            out_vec.push_back(temp_vec[0]);
            temp_vec.erase(temp_vec.begin());
            if (temp_vec.size())
            {
                if (comp_q == temp_vec[0]->get_quality())
                    std::sort(temp_vec.begin(),temp_vec.end(),sortRule2);
                temp_vec[0]->set_tracknumber(2);
                out_vec.push_back(temp_vec[0]);
                temp_vec.erase(temp_vec.begin());
            }
        }
        wipe(temp_vec);
    }
}

int getTempALCTQuality(int temp_quality)
{
    return (temp_quality > 3) ? temp_quality - 3 : 0; 
}

void wipe(std::vector<ALCTCandidate*> candvec)
{
    for (int i=0; i<candvec.size(); i++)
    {
        delete candvec[i];
    }
}

void wipe(std::vector<ALCT_ChamberHits*>cvec)
{
    for (int i=0; i<cvec.size(); i++)
    {
        delete cvec[i];
    }
}

