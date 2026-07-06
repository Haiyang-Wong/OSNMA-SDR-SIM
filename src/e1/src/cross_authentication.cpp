#include "../include/galileo-sdr.h"
#include <algorithm>
#include <climits>
#include <sstream>

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

bool isPageVisibleForPrn(int prn, int tow)
{
    if (tow_visible_prn_sets.empty()) {
        return true;
    }

    for (size_t i = 0; i < tow_visible_prn_sets.size(); i++) {
        int rangeStart = tow_visible_prn_sets[i].first;
        int rangeEnd = (i + 1 < tow_visible_prn_sets.size())
                           ? tow_visible_prn_sets[i + 1].first
                           : INT_MAX;
        if (tow >= rangeStart && tow < rangeEnd) {
            const vector<int> &visiblePrns = tow_visible_prn_sets[i].second;
            return find(visiblePrns.begin(), visiblePrns.end(), prn) != visiblePrns.end();
        }
    }

    return false;
}

vector<page_struct> getPagesForPrn(const vector<page_struct> &pages, int prn)
{
    vector<page_struct> result;
    for (size_t i = 0; i < pages.size(); i++) {
        if (pages[i].prn == prn && isPageVisibleForPrn(prn, static_cast<int>(pages[i].tow))) {
            result.push_back(pages[i]);
        }
    }
    return result;
}

vector<page_struct> getAllPagesForPrn(const vector<page_struct> &pages, int prn)
{
    vector<page_struct> result;
    for (size_t i = 0; i < pages.size(); i++) {
        if (pages[i].prn == prn) {
            result.push_back(pages[i]);
        }
    }
    return result;
}

vector<string> getCrossAuthDataADKD0FromSourcePages(int prn, vector<page_struct> pages)
{
    vector<gst> gstList = getGst(pages);
    vector<vector<string> > words(max(numSub, 1), vector<string>(5));
    int subframeIndex = 0;

    for (size_t i = 0; i < pages.size(); i++) {
        int tow = static_cast<int>(pages[i].tow);
        if (subframeIndex >= static_cast<int>(words.size())) {
            break;
        }

        string navBin = hexToBinary(pages[i].navMsg);
        if (getPageNumInSub(tow) == 10) {
            words[subframeIndex][0] = navBin.substr(8, 106) + navBin.substr(122, 14);
        }
        if (getPageNumInSub(tow) == 0) {
            words[subframeIndex][1] = navBin.substr(8, 106) + navBin.substr(122, 14);
        }
        if (getPageNumInSub(tow) == 11) {
            words[subframeIndex][2] = navBin.substr(8, 106) + navBin.substr(122, 16);
        }
        if (getPageNumInSub(tow) == 1) {
            words[subframeIndex][3] = navBin.substr(8, 106) + navBin.substr(122, 14);
        }
        if (getPageNumInSub(tow) == 12) {
            words[subframeIndex][4] = navBin.substr(8, 67);
            subframeIndex++;
        }
    }

    size_t usableSubframes = min(static_cast<size_t>(max(numSub - 1, 0)),
                                 min(words.size(), gstList.size()));
    vector<string> authData;
    for (size_t i = 0; i < usableSubframes; i++) {
        string navDataBin;
        for (size_t j = 0; j < words[i].size(); j++) {
            navDataBin += words[i][j];
        }
        if (navDataBin.size() != 549) {
            authData.push_back("");
            continue;
        }

        string navStatus = dic == 33 ? "01" : "10";
        bitset<32> gstBits((stoi(gstList[i].wn) << 20) | stoi(gstList[i].tow));
        authData.push_back(bitset<8>(prn).to_string() +
                           gstBits.to_string() +
                           hexToBinary("01") +
                           navStatus +
                           navDataBin);
    }

    return authData;
}

const vector<page_struct> &getCachedSourcePagesForPrn(int prn)
{
    static map<int, vector<page_struct> > cache;
    map<int, vector<page_struct> >::iterator it = cache.find(prn);
    if (it == cache.end()) {
        it = cache.insert(make_pair(prn, getAllPagesForPrn(nav_page_all, prn))).first;
    }
    return it->second;
}

const vector<string> &getCachedAuthDataADKD0(int prn)
{
    static map<int, vector<string> > cache;
    map<int, vector<string> >::iterator it = cache.find(prn);
    if (it == cache.end()) {
        vector<page_struct> pages = getCachedSourcePagesForPrn(prn);
        it = cache.insert(make_pair(prn, getCrossAuthDataADKD0FromSourcePages(prn, pages))).first;
    }
    return it->second;
}

const vector<string> &getCachedVisibleAuthDataADKD0(int prn)
{
    static map<int, vector<string> > cache;
    map<int, vector<string> >::iterator it = cache.find(prn);
    if (it == cache.end()) {
        vector<page_struct> pages = getCachedSourcePagesForPrn(prn);
        it = cache.insert(make_pair(prn, getAuthData(prn, pages))).first;
    }
    return it->second;
}

string makeAdkd12CacheKey(int targetPrn, int neighborPrn, int tagIndex)
{
    ostringstream key;
    key << targetPrn << ":" << neighborPrn << ":" << tagIndex;
    return key.str();
}

const vector<string> &getCachedCrossAuthDataADKD12(int targetPrn, int neighborPrn, int tagIndex);

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
    string value;
    value.push_back(kPrnMapping[(prn >> 4) & 0x0F]);
    value.push_back(kPrnMapping[prn & 0x0F]);
    return value;
}

string makeTagInfoBits(int prn_d)
{
    return hexToBinary(prnToOsnmaHex(prn_d)) + "0000" + "0001";
}

string makeTagInfoBits(int prn_d, int adkd)
{
    return hexToBinary(prnToOsnmaHex(prn_d)) + bitset<4>(adkd).to_string() + "0001";
}

string makeMacseqTagInfoBits(int prn_d)
{
    return makeTagInfoBits(prn_d);
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

bool readBits(const vector<page_struct> &pages, int index, int start, int len, string &bits)
{
    bits.clear();
    if (index < 0 || index >= static_cast<int>(pages.size())) {
        return false;
    }

    string navBin = hexToBinary(pages[index].navMsg);
    if (start < 0 || start + len > static_cast<int>(navBin.size())) {
        return false;
    }

    bits = navBin.substr(start, len);
    return true;
}

bool readCrossTagInfo(const vector<page_struct> &pages, int subframe_start, const CrossTagSlot &slot, string &infoBits)
{
    infoBits.clear();
    string firstBits;
    int infoIndex = subframe_start + slot.info_page - 1;
    if (!readBits(pages, infoIndex, slot.info_bit, slot.info_len, firstBits)) {
        return false;
    }
    infoBits += firstBits;
    if (slot.info_second_len > 0) {
        string secondBits;
        int secondInfoIndex = subframe_start + slot.info_second_page - 1;
        if (!readBits(pages, secondInfoIndex, slot.info_second_bit, slot.info_second_len, secondBits)) {
            return false;
        }
        infoBits += secondBits;
    }
    return infoBits.size() == 16;
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

bool writeCrossTagInfo(vector<page_struct> &pages, int subframe_start, const CrossTagSlot &slot, int prn_d, int adkd)
{
    string infoBits = makeTagInfoBits(prn_d, adkd);
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

bool clearCrossTagSlot(vector<page_struct> &pages, int subframe_start, const CrossTagSlot &slot)
{
    bool ok = writeBits(pages, subframe_start + slot.first_page - 1, slot.first_bit, slot.first_len, string(slot.first_len, '0')) &&
              writeBits(pages, subframe_start + slot.second_page - 1, slot.second_bit, slot.second_len, string(slot.second_len, '0')) &&
              writeBits(pages, subframe_start + slot.info_page - 1, slot.info_bit, slot.info_len, string(slot.info_len, '0'));
    if (slot.info_second_len > 0) {
        ok = ok && writeBits(pages,
                             subframe_start + slot.info_second_page - 1,
                             slot.info_second_bit,
                             slot.info_second_len,
                             string(slot.info_second_len, '0'));
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

string findKeyForFirstTow(const vector<KeyItem> &keyItems, int firstTow)
{
    for (size_t i = 0; i < keyItems.size(); i++) {
        if (keyItems[i].firstTow == firstTow) {
            return keyItems[i].key;
        }
    }
    return "";
}

vector<int> makeNeighborQueue(const vector<int> &neighbors, int pair_index)
{
    vector<int> queue7;
    if (neighbors.empty()) {
        return queue7;
    }

    int usable = min(neb_num, static_cast<int>(neighbors.size()));
    int startIndex = static_cast<int>((static_cast<long long>(pair_index) * 6) % usable);
    for (int i = 0; i < 6; i++) {
        int index = (startIndex + i) % usable;
        queue7.push_back(neighbors[index]);
    }
    return queue7;
}

string formatPrnList(const vector<int> &prns)
{
    ostringstream out;
    for (size_t i = 0; i < prns.size(); i++) {
        if (i > 0) {
            out << ",";
        }
        out << prns[i];
    }
    return out.str();
}

string formatPrnList(const vector<int> &prns, size_t limit)
{
    vector<int> limited;
    for (size_t i = 0; i < prns.size() && i < limit; i++) {
        limited.push_back(prns[i]);
    }
    return formatPrnList(limited);
}

bool hasPrn(const vector<int> &prns, int prn)
{
    return find(prns.begin(), prns.end(), prn) != prns.end();
}

void appendUniquePrn(vector<int> &prns, int prn)
{
    if (!hasPrn(prns, prn)) {
        prns.push_back(prn);
    }
}

struct NeighborSelectionSkips {
    vector<int> self;
    vector<int> duplicate;
    vector<int> osnmaEnabled;
    vector<int> disallowed;
    vector<int> missingNav;
};

bool addSourceNeighborCandidate(int targetPrn,
                                int candidatePrn,
                                int targetCount,
                                vector<int> &selectedNeighbors,
                                NeighborSelectionSkips &skips)
{
    if (targetCount > 0 && static_cast<int>(selectedNeighbors.size()) >= targetCount) {
        return true;
    }
    if (candidatePrn == targetPrn) {
        appendUniquePrn(skips.self, candidatePrn);
        return false;
    }
    if (hasPrn(selectedNeighbors, candidatePrn)) {
        appendUniquePrn(skips.duplicate, candidatePrn);
        return false;
    }
    if (!canBeCrossAuthSourcePrnd(targetPrn, candidatePrn, cross_auth_mode)) {
        if (cross_auth_mode == 3 && isOsnmaEnabledPrn(candidatePrn)) {
            appendUniquePrn(skips.osnmaEnabled, candidatePrn);
        } else {
            appendUniquePrn(skips.disallowed, candidatePrn);
        }
        return false;
    }
    if (getCachedSourcePagesForPrn(candidatePrn).empty()) {
        appendUniquePrn(skips.missingNav, candidatePrn);
        return false;
    }

    appendUniquePrn(selectedNeighbors, candidatePrn);
    return true;
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
                 const vector<CrossTagSlot> &flexSlots,
                 const vector<KeyItem> &targetKeyItems,
                 int targetPrn)
{
    (void)subframeIndex;
    int tow = static_cast<int>(targetPages[subframeStart].tow);
    int firstTow = getFirstTow(tow);
    string key = findKeyForFirstTow(targetKeyItems, firstTow + 30);
    if (key.empty()) {
        return false;
    }

    bitset<32> gstBits((targetPages[subframeStart].wn << 20) | (tow - 1));
    string authData = hexToBinary(prnToOsnmaHex(targetPrn)) + gstBits.to_string();
    for (size_t i = 0; i < flexSlots.size(); i++) {
        string infoBits;
        if (!readCrossTagInfo(targetPages, subframeStart, flexSlots[i], infoBits)) {
            return false;
        }
        authData += infoBits;
    }

    string macseqBits = hexToBinary(computeMacseq(key, authData)).substr(0, 12);
    int page2Index = subframeStart + 1;
    return writeBits(targetPages, page2Index, 154, 12, macseqBits);
}

bool hasValidFlexPrn(const vector<int> &flexPrns)
{
    for (size_t i = 0; i < flexPrns.size(); i++) {
        if (flexPrns[i] != 0) {
            return true;
        }
    }
    return false;
}

int extractSelfAdkd0AuthTow(const string &authData)
{
    if (authData.size() < 40) {
        return -1;
    }
    return static_cast<int>(strtol(authData.substr(20, 20).c_str(), NULL, 2));
}

int findAuthDataIndexByTow(const vector<string> &authDataList, int tagGstTow)
{
    for (size_t i = 0; i < authDataList.size(); i++) {
        if (extractSelfAdkd0AuthTow(authDataList[i]) == tagGstTow) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

int extractAdkd12AuthTow(const string &authData)
{
    if (authData.size() < 48) {
        return -1;
    }
    return static_cast<int>(strtol(authData.substr(28, 20).c_str(), NULL, 2));
}

int findAdkd12AuthDataIndexByTow(const vector<string> &authDataList, int tagGstTow)
{
    for (size_t i = 0; i < authDataList.size(); i++) {
        if (extractAdkd12AuthTow(authDataList[i]) == tagGstTow) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

bool writeSlotForNeighbor(vector<page_struct> &targetPages,
                          int subframeStart,
                          int subframeIndex,
                          const CrossTagSlot &slot,
                          int targetPrn,
                          int neighborPrn,
                          const vector<KeyItem> &targetKeyItems)
{
    int authIndex = subframeIndex - 1;
    if (authIndex < 0) {
        return false;
    }
    int firstTow = getFirstTow(static_cast<int>(targetPages[subframeStart].tow));
    int tagGstTow = firstTow - 1;
    string key = findKeyForFirstTow(targetKeyItems, firstTow + 30);
    if (key.empty()) {
        return false;
    }
    const vector<page_struct> &neighborPages = getCachedSourcePagesForPrn(neighborPrn);
    if (neighborPages.empty()) {
        return false;
    }

    int writeTow = static_cast<int>(targetPages[subframeStart].tow);
    bool prndVisibleAtWriteTow = isPageVisibleForPrn(neighborPrn, writeTow);
    const vector<string> &neighborAuthData = prndVisibleAtWriteTow
                                                ? getCachedVisibleAuthDataADKD0(neighborPrn)
                                                : getCachedAuthDataADKD0(neighborPrn);
    int authDataIndex = findAuthDataIndexByTow(neighborAuthData, tagGstTow);
    if (authDataIndex < 0) {
        return false;
    }

    const string &sourceAuthData = neighborAuthData[authDataIndex];
    string authData = makeCrossAuthData(targetPrn, slot.tag_index, sourceAuthData);
    string tag = computeSingleCrossTag(key, authData);
    bool ok = writeCrossTag(targetPages, subframeStart, slot, tag) &&
              writeCrossTagInfo(targetPages, subframeStart, slot, neighborPrn);
    if (ok) {
        recordGeneratedOsnmaTag(writeTow,
                                tagGstTow,
                                targetPrn,
                                neighborPrn,
                                0,
                                targetPrn != neighborPrn,
                                "cross_adkd0",
                                slot.tag_index);
    }
    return ok;
}

vector<string> getCrossAuthDataADKD12(int targetPrn,
                                      int neighborPrn,
                                      int tagIndex,
                                      const vector<page_struct> &neighborPages)
{
    vector<string> authData;
    vector<gst> gstList = getGstADKD12(neighborPages);
    vector<string> navDataBySubframe;
    vector<string> words(5);

    for (size_t i = 0; i < neighborPages.size(); i++) {
        int tow = static_cast<int>(neighborPages[i].tow);
        string navBin = hexToBinary(neighborPages[i].navMsg);
        if (getPageNumInSub(tow) == 10) {
            words[0] = navBin.substr(8, 106) + navBin.substr(122, 14);
        }
        if (getPageNumInSub(tow) == 0) {
            words[1] = navBin.substr(8, 106) + navBin.substr(122, 14);
        }
        if (getPageNumInSub(tow) == 11) {
            words[2] = navBin.substr(8, 106) + navBin.substr(122, 16);
        }
        if (getPageNumInSub(tow) == 1) {
            words[3] = navBin.substr(8, 106) + navBin.substr(122, 14);
        }
        if (getPageNumInSub(tow) == 12) {
            words[4] = navBin.substr(8, 67);
            string navDataBin;
            for (size_t j = 0; j < words.size(); j++) {
                navDataBin += words[j];
            }
            if (navDataBin.size() == 549) {
                navDataBySubframe.push_back(navDataBin);
            }
            words = vector<string>(5);
        }
    }
    if (navDataBySubframe.empty()) {
        return authData;
    }

    string prnD = hexToBinary(prnToOsnmaHex(neighborPrn));
    string prnA = hexToBinary(prnToOsnmaHex(targetPrn));
    string counter = bitset<8>(tagIndex).to_string();
    string navStatus = dic == 33 ? "01" : "10";

    for (size_t i = 0; i < gstList.size(); i++) {
        size_t navDataIndex = min(i, navDataBySubframe.size() - 1);
        bitset<32> gstBits((stoi(gstList[i].wn) << 20) | stoi(gstList[i].tow));
        authData.push_back(prnD + prnA + gstBits.to_string() + counter + navStatus + navDataBySubframe[navDataIndex]);
    }
    return authData;
}

const vector<string> &getCachedCrossAuthDataADKD12(int targetPrn, int neighborPrn, int tagIndex)
{
    static map<string, vector<string> > cache;
    string key = makeAdkd12CacheKey(targetPrn, neighborPrn, tagIndex);
    map<string, vector<string> >::iterator it = cache.find(key);
    if (it == cache.end()) {
        const vector<page_struct> &pages = getCachedSourcePagesForPrn(neighborPrn);
        it = cache.insert(make_pair(key, getCrossAuthDataADKD12(targetPrn, neighborPrn, tagIndex, pages))).first;
    }
    return it->second;
}

bool writeSlotForNeighborADKD12(vector<page_struct> &targetPages,
		                                int subframeStart,
		                                int subframeIndex,
	                                const CrossTagSlot &slot,
	                                int targetPrn,
	                                int neighborPrn,
	                                const vector<KeyItem> &targetKeyItems)
	{
	    if (subframeIndex <= 0) {
	        return false;
	    }
	    int firstTow = getFirstTow(static_cast<int>(targetPages[subframeStart].tow));
	    int tagGstTow = firstTow - 1;
	    string key = findKeyForFirstTow(targetKeyItems, firstTow + 330);
	    if (key.empty()) {
	        return false;
	    }
    const vector<page_struct> &neighborPages = getCachedSourcePagesForPrn(neighborPrn);
	    if (neighborPages.empty()) {
	        return false;
	    }
	
	    const vector<string> &neighborAuthData = getCachedCrossAuthDataADKD12(targetPrn, neighborPrn, slot.tag_index);
	    int authIndex = findAdkd12AuthDataIndexByTow(neighborAuthData, tagGstTow);
	    if (authIndex < 0) {
	        return false;
		    }
	
		    string tag = computeSingleCrossTag(key, neighborAuthData[authIndex]);
	    bool ok = writeCrossTag(targetPages, subframeStart, slot, tag) &&
	              writeCrossTagInfo(targetPages, subframeStart, slot, neighborPrn, 12);
		    if (ok) {
		        int writeTow = static_cast<int>(targetPages[subframeStart].tow);
		        recordGeneratedOsnmaTag(writeTow,
		                                tagGstTow,
		                                targetPrn,
		                                neighborPrn,
	                                12,
	                                targetPrn != neighborPrn,
	                                "cross_adkd12",
	                                slot.tag_index);
	    }
	    return ok;
	}

bool writeSlotForAnyNeighbor(vector<page_struct> &targetPages,
		                             int subframeStart,
		                             int subframeIndex,
		                             const CrossTagSlot &slot,
		                             int targetPrn,
		                             const vector<int> &preferredQueue,
		                             int preferredOffset,
		                             const vector<int> &fallbackNeighbors,
		                             const vector<KeyItem> &targetKeyItems,
		                             vector<int> &usedPrns,
		                             int &writtenPrn)
{
    vector<int> candidates;
    for (int i = 0; i < static_cast<int>(preferredQueue.size()); i++) {
        int index = (preferredOffset + i) % static_cast<int>(preferredQueue.size());
        appendUniquePrn(candidates, preferredQueue[index]);
    }
    for (size_t i = 0; i < fallbackNeighbors.size(); i++) {
        appendUniquePrn(candidates, fallbackNeighbors[i]);
    }

    for (size_t i = 0; i < candidates.size(); i++) {
        int prnd = candidates[i];
        if (fallbackNeighbors.size() >= 3 && hasPrn(usedPrns, prnd)) {
            continue;
        }
        if (writeSlotForNeighbor(targetPages, subframeStart, subframeIndex, slot, targetPrn, prnd, targetKeyItems)) {
            appendUniquePrn(usedPrns, prnd);
            writtenPrn = prnd;
            return true;
        }
    }

    writtenPrn = 0;
    return false;
}

bool writeSlotForAnyNeighborADKD12(vector<page_struct> &targetPages,
		                                   int subframeStart,
		                                   int subframeIndex,
	                                   const CrossTagSlot &slot,
	                                   int targetPrn,
	                                   const vector<int> &preferredQueue,
	                                   int preferredOffset,
	                                   const vector<int> &fallbackNeighbors,
	                                   const vector<KeyItem> &targetKeyItems,
	                                   vector<int> &usedPrns,
	                                   int &writtenPrn)
{
    vector<int> candidates;
    for (int i = 0; i < static_cast<int>(preferredQueue.size()); i++) {
        int index = (preferredOffset + i) % static_cast<int>(preferredQueue.size());
        appendUniquePrn(candidates, preferredQueue[index]);
    }
    for (size_t i = 0; i < fallbackNeighbors.size(); i++) {
        appendUniquePrn(candidates, fallbackNeighbors[i]);
    }

    for (size_t i = 0; i < candidates.size(); i++) {
        int prnd = candidates[i];
        if (hasPrn(usedPrns, prnd)) {
            continue;
        }
        if (writeSlotForNeighborADKD12(targetPages, subframeStart, subframeIndex, slot, targetPrn, prnd, targetKeyItems)) {
            usedPrns.push_back(prnd);
            writtenPrn = prnd;
            return true;
        }
    }

    writtenPrn = 0;
    return false;
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

    int targetNeighborCount = neb_num > 0 ? neb_num : static_cast<int>(neighbors.size());

    vector<int> selectedNeighbors;
    NeighborSelectionSkips skips;
    for (size_t i = 0; i < neighbors.size(); i++) {
        addSourceNeighborCandidate(target_prn, neighbors[i], targetNeighborCount, selectedNeighbors, skips);
        if (static_cast<int>(selectedNeighbors.size()) >= targetNeighborCount) {
            break;
        }
    }

    if (cross_auth_mode == 3) {
        for (size_t i = 0; i < no_osnma_prns.size(); i++) {
            addSourceNeighborCandidate(target_prn, no_osnma_prns[i], targetNeighborCount, selectedNeighbors, skips);
            if (static_cast<int>(selectedNeighbors.size()) >= targetNeighborCount) {
                break;
            }
        }
    }

    cout << "[CrossAuth] E" << target_prn
         << " authenticated with PRNd candidates: [" << formatPrnList(selectedNeighbors) << "]" << endl;

    bool hasCrossNeighbors = !selectedNeighbors.empty();

    vector<page_struct> targetPages = getPagesForPrn(nav_page, target_prn);
    if (targetPages.empty()) {
        return;
    }
    normalizePagesToHex(targetPages);

    vector<KeyItem> targetKeyItems = getPageKeyItems(targetPages);
    if (targetKeyItems.empty()) {
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
        vector<int> queue7;
        if (hasCrossNeighbors) {
            queue7 = makeNeighborQueue(selectedNeighbors, pairIndex);
            if (queue7.size() != 6) {
                continue;
            }
        }

        if (is60page(firstTow)) {
            bool cleared = clearCrossTagSlot(targetPages, start, tag2) &&
                           clearCrossTagSlot(targetPages, start, tag4) &&
                           clearCrossTagSlot(targetPages, start, tag6);
            changed |= cleared;
            bool tag2Written = false;
            bool tag4Written = false;
            bool tag6Written = false;
            int tag2Prn = 0;
            int tag4Prn = 0;
            vector<int> usedPrns;
            if (hasCrossNeighbors && subframeIndex > 0) {
                tag2Written = writeSlotForAnyNeighbor(targetPages, start, subframeIndex, tag2, target_prn,
	                                                      queue7, 0, selectedNeighbors, targetKeyItems, usedPrns, tag2Prn);
	                tag4Written = writeSlotForAnyNeighbor(targetPages, start, subframeIndex, tag4, target_prn,
	                                                      queue7, 1, selectedNeighbors, targetKeyItems, usedPrns, tag4Prn);
                changed |= tag2Written;
                changed |= tag4Written;
                int tag6Prn = 0;
                tag6Written = writeSlotForAnyNeighbor(targetPages, start, subframeIndex, tag6, target_prn,
	                                                      queue7, 2, selectedNeighbors, targetKeyItems, usedPrns, tag6Prn);
                changed |= tag6Written;
            }
	            if (tag2Written && tag4Written && tag6Written) {
	                vector<CrossTagSlot> flexSlots;
	                flexSlots.push_back(tag2);
	                flexSlots.push_back(tag4);
	                changed |= writeMacseq(targetPages, start, subframeIndex, flexSlots, targetKeyItems, target_prn);
	            }
        } else {
            bool cleared = clearCrossTagSlot(targetPages, start, tag2) &&
                           clearCrossTagSlot(targetPages, start, tag3) &&
                           clearCrossTagSlot(targetPages, start, tag5) &&
                           clearCrossTagSlot(targetPages, start, tag6);
            changed |= cleared;
            bool tag2Written = false;
            int tag2Prn = 0;
            vector<int> usedPrns;
            if (hasCrossNeighbors && subframeIndex > 0) {
                tag2Written = writeSlotForAnyNeighbor(targetPages, start, subframeIndex, tag2, target_prn,
	                                                      queue7, 3, selectedNeighbors, targetKeyItems, usedPrns, tag2Prn);
                changed |= tag2Written;
                int tag3Prn = 0;
                int tag5Prn = 0;
                int tag6Prn = 0;
                changed |= writeSlotForAnyNeighbor(targetPages, start, subframeIndex, tag3, target_prn,
	                                                   queue7, 4, selectedNeighbors, targetKeyItems, usedPrns, tag3Prn);
                changed |= writeSlotForAnyNeighbor(targetPages, start, subframeIndex, tag5, target_prn,
	                                                   queue7, 5, selectedNeighbors, targetKeyItems, usedPrns, tag5Prn);
                if (!selectedNeighbors.empty()) {
                    changed |= writeSlotForNeighborADKD12(targetPages, start, subframeIndex, tag6, target_prn,
                                                          selectedNeighbors.front(), targetKeyItems);
                    tag6Prn = selectedNeighbors.front();
                }
            }
	            if (tag2Written) {
	                vector<CrossTagSlot> flexSlots;
	                flexSlots.push_back(tag2);
	                changed |= writeMacseq(targetPages, start, subframeIndex, flexSlots, targetKeyItems, target_prn);
	            }
        }
    }

    if (changed) {
        fillCrc(targetPages);
        singleIn(target_prn, targetPages);
    }
}
