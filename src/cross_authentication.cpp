#include "../include/galileo-sdr.h"

namespace {

const string kPrnMapping = "0123456789ABCDEF";

struct CrossTagSlot {
    int tag_index;
    int first_page;
    int first_bit;
    int first_len;
    int second_page;
    int second_bit;
    int second_len;
    int info_page;
    int info_bit;
    int info_len;
    int info_second_page;
    int info_second_bit;
    int info_second_len;
};

vector<page_struct> getPagesForPrn(const vector<page_struct> &pages, int prn)
{
    vector<page_struct> result;
    for (size_t i = 0; i < pages.size(); i++) {
        if (pages[i].prn == prn) {
            result.push_back(pages[i]);
        }
    }
    return result;
}

bool isBinaryPageString(const string &value)
{
    if (value.size() != 240) {
        return false;
    }
    for (size_t i = 0; i < value.size(); i++) {
        if (value[i] != '0' && value[i] != '1') {
            return false;
        }
    }
    return true;
}

void normalizePagesToHex(vector<page_struct> &pages)
{
    for (size_t i = 0; i < pages.size(); i++) {
        if (isBinaryPageString(pages[i].navMsg)) {
            pages[i].navMsg = binaryToHexBitset(pages[i].navMsg);
        }
    }
}

string prnToOsnmaHex(int prn)
{
    if (prn < 16) {
        return string("0") + kPrnMapping[prn];
    }
    if (prn < 30) {
        return to_string(prn - 6);
    }
    return to_string(prn - 12);
}

string makeTagInfoBits(int prn_d)
{
    return hexToBinary(prnToOsnmaHex(prn_d)) + "0000" + "0001";
}

string makeMacseqTagInfoBits(int prn_d)
{
    //if (prn_d == 0) {
     //   return string(16, '0');
    //}
   // return makeTagInfoBits(prn_d);
}

string makeCrossAuthData(int targetPrn, int tagIndex, const string &neighborAuthData)
{
    if (neighborAuthData.size() < 48) {
        return neighborAuthData;
    }
    return neighborAuthData.substr(0, 8) +
           hexToBinary(prnToOsnmaHex(targetPrn)) +
           neighborAuthData.substr(8, 32) +
           bitset<8>(tagIndex).to_string() +
           neighborAuthData.substr(48);
}

bool writeBits(vector<page_struct> &pages, int index, int start, int len, const string &bits)
{
    if (index < 0 || index >= static_cast<int>(pages.size()) || static_cast<int>(bits.size()) != len) {
        return false;
    }

    string navBin = hexToBinary(pages[index].navMsg);
    if (start < 0 || start + len > static_cast<int>(navBin.size())) {
        return false;
    }

    navBin.replace(start, len, bits);
    pages[index].navMsg = binaryToHexBitset(navBin);
    return true;
}

bool writeCrossTag(vector<page_struct> &pages, int subframe_start, const CrossTagSlot &slot, const string &tag_hex)
{
    string tagBits = hexToBinary(tag_hex.substr(0, 10));
    if (tagBits.size() != 40) {
        return false;
    }

    int firstIndex = subframe_start + slot.first_page - 1;
    int secondIndex = subframe_start + slot.second_page - 1;
    string firstBits = tagBits.substr(0, slot.first_len);
    string secondBits = tagBits.substr(slot.first_len, slot.second_len);

    return writeBits(pages, firstIndex, slot.first_bit, slot.first_len, firstBits) &&
           writeBits(pages, secondIndex, slot.second_bit, slot.second_len, secondBits);
}

bool writeCrossTagInfo(vector<page_struct> &pages, int subframe_start, const CrossTagSlot &slot, int prn_d)
{
    string infoBits = makeTagInfoBits(prn_d);
    int infoIndex = subframe_start + slot.info_page - 1;
    bool ok = writeBits(pages, infoIndex, slot.info_bit, slot.info_len, infoBits.substr(0, slot.info_len));
    if (slot.info_second_len > 0) {
        int secondInfoIndex = subframe_start + slot.info_second_page - 1;
        ok = ok && writeBits(pages,
                             secondInfoIndex,
                             slot.info_second_bit,
                             slot.info_second_len,
                             infoBits.substr(slot.info_len, slot.info_second_len));
    }
    return ok;
}

string computeSingleCrossTag(const string &key_hex, const string &auth_data)
{
    vector<uint8_t> keyBytes = hexstr_to_bytes(key_hex);
    vector<uint8_t> dataBytes = bitstring_to_bytes(auth_data);
    vector<uint8_t> mac = hmac_sha256(keyBytes, dataBytes);
    return macToString(mac).substr(0, 10);
}

string computeMacseq(const string &key_hex, const string &auth_data)
{
    vector<uint8_t> keyBytes = hexstr_to_bytes(key_hex);
    vector<uint8_t> dataBytes = bitstring_to_bytes(auth_data);
    vector<uint8_t> mac = hmac_sha256(keyBytes, dataBytes);
    return macToString(mac).substr(0, 4);
}

vector<int> makeNeighborQueue(const vector<int> &neighbors, int pair_index)
{
    vector<int> queue7;
    if (neighbors.empty()) {
        return queue7;
    }

    int usable = min(neb_num, static_cast<int>(neighbors.size()));
    for (int i = 0; i < 7; i++) {
        int index = (pair_index + i) % usable;
        queue7.push_back(neighbors[index]);
    }
    return queue7;
}

bool pageMatchesSubframeOffset(const vector<page_struct> &pages, int index, int expected_page_num)
{
    if (index < 0 || index >= static_cast<int>(pages.size())) {
        return false;
    }
    int tow = static_cast<int>(pages[index].tow);
    return getPageNumInSub(tow) == expected_page_num;
}

bool writeMacseq(vector<page_struct> &targetPages,
                 int subframeStart,
                 int subframeIndex,
                 const vector<int> &flexPrns,
                 const vector<string> &targetKeys,
                 int targetPrn)
{
    int keyIndex = subframeIndex + 1;
    if (keyIndex < 0 || keyIndex >= static_cast<int>(targetKeys.size())) {
        return false;
    }

    int tow = static_cast<int>(targetPages[subframeStart].tow);
    bitset<32> gstBits((targetPages[subframeStart].wn << 20) | (tow - 1));
    string authData = hexToBinary(prnToOsnmaHex(targetPrn)) + gstBits.to_string();
    for (size_t i = 0; i < flexPrns.size(); i++) {
        authData += makeTagInfoBits(flexPrns[i]);

    }

    string macseqBits = hexToBinary(computeMacseq(targetKeys[keyIndex], authData)).substr(0, 12);
    int page2Index = subframeStart + 1;
    return writeBits(targetPages, page2Index, 154, 12, macseqBits);
}

bool writeSlotForNeighbor(vector<page_struct> &targetPages,
                          int subframeStart,
                          int subframeIndex,
                          const CrossTagSlot &slot,
                          int targetPrn,
                          int neighborPrn,
                          const vector<string> &targetKeys)
{
    int authIndex = subframeIndex - 1;
    int keyIndex = subframeIndex + 1;
    if (authIndex < 0 || keyIndex < 0 || keyIndex >= static_cast<int>(targetKeys.size())) {
        return false;
    }

    vector<page_struct> neighborPages = getPagesForPrn(nav_page_all, neighborPrn);
    if (neighborPages.empty()) {
        return false;
    }

    vector<string> neighborAuthData = getAuthData(neighborPrn, neighborPages);
    if (authIndex >= static_cast<int>(neighborAuthData.size())) {
        return false;
    }

    string authData = makeCrossAuthData(targetPrn, slot.tag_index, neighborAuthData[authIndex]);
    string tag = computeSingleCrossTag(targetKeys[keyIndex], authData);
    return writeCrossTag(targetPages, subframeStart, slot, tag) &&
           writeCrossTagInfo(targetPages, subframeStart, slot, neighborPrn);
}

} // namespace

void runCrossAuthenticationForVisiblePrn(int visible_index, int target_prn)
{
    if (visible_index < 0 || visible_index >= static_cast<int>(nearest_prn_by_visible.size())) {
        return;
    }

    const vector<int> &neighbors = nearest_prn_by_visible[visible_index];
    if (neighbors.empty()) {
        return;
    }

    vector<page_struct> targetPages = getPagesForPrn(nav_page, target_prn);
    if (targetPages.empty()) {
        return;
    }
    normalizePagesToHex(targetPages);

    vector<string> targetKeys = getKey(targetPages);
    if (targetKeys.empty()) {
        return;
    }

    const CrossTagSlot tag2 = {2, 2, 170, 8, 3, 146, 32, 4, 146, 16, 0, 0, 0};
    const CrossTagSlot tag3 = {3, 4, 162, 16, 5, 146, 24, 5, 170, 8, 6, 146, 8};
    const CrossTagSlot tag4 = {4, 6, 154, 24, 7, 146, 16, 7, 162, 16, 0, 0, 0};
    const CrossTagSlot tag5 = {5, 8, 146, 32, 9, 146, 8, 9, 154, 16, 0, 0, 0};
    const CrossTagSlot tag6 = {6, 9, 170, 8, 10, 146, 32, 11, 146, 16, 0, 0, 0};

    bool changed = false;
    for (int start = 0; start < static_cast<int>(targetPages.size()); start++) {
        int tow = static_cast<int>(targetPages[start].tow);
        if (getPageNumInSub(tow) != 0) {
            continue;
        }

        int firstTow = getFirstTow(tow);
        int subframeIndex = (firstTow - static_cast<int>(start_time)) / 30;

        if (!pageMatchesSubframeOffset(targetPages, start + 10, 10)) {
            continue;
        }

        int pairIndex = (firstTow - static_cast<int>(start_time)) / 60;
        vector<int> queue7 = makeNeighborQueue(neighbors, pairIndex);
        if (queue7.size() != 7) {
            continue;
        }

        if (is60page(firstTow)) {
            if (subframeIndex > 0) {
                changed |= writeSlotForNeighbor(targetPages, start, subframeIndex, tag2, target_prn, queue7[0], targetKeys);
                changed |= writeSlotForNeighbor(targetPages, start, subframeIndex, tag4, target_prn, queue7[1], targetKeys);
                changed |= writeSlotForNeighbor(targetPages, start, subframeIndex, tag6, target_prn, queue7[2], targetKeys);
            }
            vector<int> flexPrns;
            if (subframeIndex > 0) {
                flexPrns.push_back(queue7[0]);
                flexPrns.push_back(queue7[1]);
            } else {
                flexPrns.push_back(0);
                flexPrns.push_back(0);
            }
            changed |= writeMacseq(targetPages, start, subframeIndex, flexPrns, targetKeys, target_prn);
        } else {
            if (subframeIndex > 0) {
                changed |= writeSlotForNeighbor(targetPages, start, subframeIndex, tag2, target_prn, queue7[3], targetKeys);
                changed |= writeSlotForNeighbor(targetPages, start, subframeIndex, tag3, target_prn, queue7[4], targetKeys);
                changed |= writeSlotForNeighbor(targetPages, start, subframeIndex, tag5, target_prn, queue7[5], targetKeys);
            }
            vector<int> flexPrns;
            if (subframeIndex > 0) {
                flexPrns.push_back(queue7[3]);
            } else {
                flexPrns.push_back(0);
            }
            changed |= writeMacseq(targetPages, start, subframeIndex, flexPrns, targetKeys, target_prn);
            // queue7[6] reserves the second-subframe Tag6/ADKD12 slot and advances the next pair.
        }
    }

    if (changed) {
        fillCrc(targetPages);
        singleIn(target_prn, targetPages);
    }
}
