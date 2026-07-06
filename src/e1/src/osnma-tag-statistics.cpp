#include "../include/galileo-sdr.h"
#include <algorithm>
#include <climits>
#include <iomanip>

vector<OsnmaTagLogRecord> osnma_tag_log_records;

namespace {

struct TagCounterStats
{
    int total = 0;
    int self = 0;
    int cross = 0;
};

struct SatelliteTagStats
{
    int total = 0;
    int self = 0;
    int cross = 0;
    map<int, int> by_adkd;
    set<int> authenticated_prnds;
    bool has_time = false;
    int first_tow = 0;
    int last_tow = 0;
};

bool isPrnVisibleAtTow(int prn, int tow)
{
    if (tow_visible_prn_sets.empty())
    {
        return true;
    }

    for (size_t i = 0; i < tow_visible_prn_sets.size(); i++)
    {
        int rangeStart = tow_visible_prn_sets[i].first;
        int rangeEnd = (i + 1 < tow_visible_prn_sets.size())
                           ? tow_visible_prn_sets[i + 1].first
                           : INT_MAX;
        if (tow >= rangeStart && tow < rangeEnd)
        {
            const vector<int> &visiblePrns = tow_visible_prn_sets[i].second;
            return find(visiblePrns.begin(), visiblePrns.end(), prn) != visiblePrns.end();
        }
    }

    return false;
}

string formatGeneratedTagUtc(int tow)
{
    galtime_t tagTime;
    tagTime.week = cur_wn + 1024;
    tagTime.sec = static_cast<double>(tow - GPS_LEAPSECONDS);
    while (tagTime.sec < 0)
    {
        tagTime.week--;
        tagTime.sec += SECONDS_IN_WEEK;
    }
    while (tagTime.sec >= SECONDS_IN_WEEK)
    {
        tagTime.week++;
        tagTime.sec -= SECONDS_IN_WEEK;
    }

    datetime_t date;
    gal2date(&tagTime, &date);
    int seconds = static_cast<int>(date.sec + 0.5);
    if (seconds >= 60)
    {
        seconds = 59;
    }

    ostringstream out;
    out << setfill('0')
        << setw(4) << date.y << "-"
        << setw(2) << date.m << "-"
        << setw(2) << date.d << "T"
        << setw(2) << date.hh << ":"
        << setw(2) << date.mm << ":"
        << setw(2) << seconds << "+00:00";
    return out.str();
}

void updateTagCounter(TagCounterStats &stats, bool is_cross)
{
    stats.total++;
    if (is_cross)
    {
        stats.cross++;
    }
    else
    {
        stats.self++;
    }
}

} // namespace

void recordGeneratedOsnmaTag(int write_tow,
                             int auth_tow,
                             int prn_a,
                             int prn_d,
                             int adkd,
                             bool is_cross,
                             const string &source,
                             int tag_slot)
{
    OsnmaTagLogRecord record;
    record.write_tow = write_tow;
    record.auth_tow = auth_tow;
    record.prn_a = prn_a;
    record.prn_d = prn_d;
    record.adkd = adkd;
    record.is_cross = is_cross;
    record.source = source;
    record.tag_slot = tag_slot;
    record.auth_type = is_cross ? "cross" : "self";
    record.prna_visible = isPrnVisibleAtTow(prn_a, write_tow);
    record.prnd_visible = isPrnVisibleAtTow(prn_d, write_tow);
    record.prnd_has_osnma = isOsnmaEnabledPrn(prn_d);
    record.prnd_invisible_cross_auth = is_cross && !record.prnd_visible;
    osnma_tag_log_records.push_back(record);
}

void writeGeneratedOsnmaTagLog(const string &output_path)
{
    ofstream out(output_path.c_str());
    if (!out.is_open())
    {
        cerr << "[WARN] Failed to write OSNMA tag log to " << output_path << endl;
        return;
    }

    int invisiblePrndCrossAuthCount = 0;

    out << "write_tow,auth_tow,prna,prnd,adkd,is_cross,source,tag_slot,slot,auth_type,prna_visible,prnd_visible,prnd_has_osnma,prnd_invisible_cross_auth\n";
    for (size_t i = 0; i < osnma_tag_log_records.size(); i++)
    {
        const OsnmaTagLogRecord &record = osnma_tag_log_records[i];
        if (record.prnd_invisible_cross_auth)
        {
            invisiblePrndCrossAuthCount++;
        }
        out << record.write_tow << ","
            << record.auth_tow << ","
            << record.prn_a << ","
            << record.prn_d << ","
            << record.adkd << ","
            << (record.is_cross ? 1 : 0) << ","
            << record.source << ","
            << record.tag_slot << ","
            << record.tag_slot << ","
            << record.auth_type << ","
            << (record.prna_visible ? 1 : 0) << ","
            << (record.prnd_visible ? 1 : 0) << ","
            << (record.prnd_has_osnma ? 1 : 0) << ","
            << (record.prnd_invisible_cross_auth ? 1 : 0) << "\n";
    }

    cout << "[OSNMA] Tag log written: " << osnma_tag_log_records.size()
         << " records -> " << output_path << endl;
    cout << "[OSNMA] Invisible PRNd cross-auth tags: "
         << invisiblePrndCrossAuthCount << endl;
}

void writeGeneratedOsnmaTagStatisticsTxt(const string &output_path)
{
    ofstream out(output_path.c_str());
    if (!out.is_open())
    {
        cerr << "[WARN] Failed to write OSNMA tag statistics to " << output_path << endl;
        return;
    }

    TagCounterStats overall;
    map<int, TagCounterStats> byAdkd;
    map<int, SatelliteTagStats> byPrna;

    for (size_t i = 0; i < osnma_tag_log_records.size(); i++)
    {
        const OsnmaTagLogRecord &record = osnma_tag_log_records[i];
        bool isCross = record.is_cross || record.prn_a != record.prn_d;
        int statTow = record.auth_tow > 0 ? record.auth_tow : record.write_tow;

        updateTagCounter(overall, isCross);
        updateTagCounter(byAdkd[record.adkd], isCross);

        SatelliteTagStats &satStats = byPrna[record.prn_a];
        satStats.total++;
        if (isCross)
        {
            satStats.cross++;
        }
        else
        {
            satStats.self++;
        }
        satStats.by_adkd[record.adkd]++;
        satStats.authenticated_prnds.insert(record.prn_d);
        if (!satStats.has_time)
        {
            satStats.has_time = true;
            satStats.first_tow = statTow;
            satStats.last_tow = statTow;
        }
        else
        {
            satStats.first_tow = min(satStats.first_tow, statTow);
            satStats.last_tow = max(satStats.last_tow, statTow);
        }
    }

    out << "===== OSNMA Tag Statistics =====\n";
    out << "Total parsed tags: " << overall.total << "\n";
    out << "SUCCESS tags: " << overall.total << "\n";
    out << "FAILURE tags: 0\n";
    out << "UNKNOWN tags: 0\n";
    out << "Skipped OSNMA lines: 0\n\n";

    out << "Self-authenticated tags: " << overall.self << "\n";
    out << "Cross-authenticated tags: " << overall.cross << "\n\n";

    out << "ADKD statistics:\n";
    const int adkdOrder[] = {0, 4, 12};
    for (size_t i = 0; i < sizeof(adkdOrder) / sizeof(adkdOrder[0]); i++)
    {
        int adkd = adkdOrder[i];
        const TagCounterStats &stats = byAdkd[adkd];
        out << "  ADKD=" << adkd
            << ": total=" << stats.total
            << ", success=" << stats.total
            << ", failure=0"
            << ", self=" << stats.self
            << ", cross=" << stats.cross << "\n";
    }

    out << "\nSatellite statistics by PRNa:\n";
    for (map<int, SatelliteTagStats>::const_iterator it = byPrna.begin(); it != byPrna.end(); ++it)
    {
        int prna = it->first;
        const SatelliteTagStats &stats = it->second;
        int adkd0 = stats.by_adkd.count(0) ? stats.by_adkd.find(0)->second : 0;
        int adkd4 = stats.by_adkd.count(4) ? stats.by_adkd.find(4)->second : 0;
        int adkd12 = stats.by_adkd.count(12) ? stats.by_adkd.find(12)->second : 0;

        out << "  E" << setfill('0') << setw(2) << prna << setfill(' ')
            << ": total=" << stats.total
            << ", self=" << stats.self
            << ", cross=" << stats.cross
            << ", ADKD0=" << adkd0
            << ", ADKD4=" << adkd4
            << ", ADKD12=" << adkd12
            << ", authenticated_PRNd_count=" << stats.authenticated_prnds.size()
            << "\n";
        if (stats.has_time)
        {
            out << "    existence: first=" << formatGeneratedTagUtc(stats.first_tow)
                << ", last=" << formatGeneratedTagUtc(stats.last_tow)
                << ", first_TOW=" << stats.first_tow
                << ", last_TOW=" << stats.last_tow
                << ", \n";
        }
    }

    cout << "[OSNMA] Tag statistics written: "
         << output_path << endl;
}
