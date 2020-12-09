#include <sdsl/suffix_arrays.hpp>
#include <string>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <vector>

using namespace sdsl;
using namespace std;

template <class t_csa, class t_rac, class t_pat_iter>
typename t_csa::size_type count_one_error_case(const t_csa &csa, t_pat_iter begin, t_pat_iter end, size_t m, bool include_middle, bool case_a, t_rac &locations, bool locate)
{
    typename t_csa::size_type x = (m + 1) / 2;
    if (!include_middle)
        x--;

    if (end - begin > (typename std::iterator_traits<t_pat_iter>::difference_type)csa.size())
        return 0;

    typename t_csa::char_type curr_char;
    typename t_csa::size_type left_res = 0, right_res = 0, left_err_res = 0, right_err_res = 0, result = 0, i = 0, occs = 0;
    size_t locations_size;
    //find SA of right half P[x...m]
    backward_search(csa, 0, csa.size() - 1, begin + x, end, left_res, right_res);

    if (left_res <= right_res)
    {
        //for each curr_char at index i in P[0...x-1]
        for (i = x; i > 0; i--)
        {
            curr_char = (typename t_csa::char_type) * (begin + i - 1);
            //check existence of P[0...i-1]<<j<<P[i+1...m] s.t j!=i, j in alphabet
            for (size_t j = 1; j < csa.sigma; j++)
            {
                if (csa.char2comp[curr_char] != j)
                {
                    backward_search(csa, left_res, right_res, csa.comp2char[j], left_err_res, right_err_res);
                    occs = backward_search(csa, left_err_res, right_err_res, begin, begin + i - 1, left_err_res, right_err_res);
                    result += occs;

                    if (locate && occs > 0)
                    {
                        locations_size = locations.size();
                        locations.resize(locations_size + occs);
                        for (typename t_csa::size_type k = 0; k < occs; k++)
                        {
                            if (case_a)
                                locations[locations_size + k] = csa[left_err_res + k];
                            else
                                locations[locations_size + k] = csa.size() - csa[left_err_res + k] - m;
                        }
                    }
                }
            }
            //return char at j index in P to the original one
            backward_search(csa, left_res, right_res, curr_char, left_res, right_res);
        }
    }
    return result;
}

template <class t_csa, class t_rac>
typename t_csa::size_type
handle_one_error(
    const t_csa &csa,
    const t_csa &rev_csa,
    string query,
    string rev_query,
    t_rac &locations,
    bool locate)
{
    size_t occs = 0;
    occs = count_one_error_case(csa, query.begin(), query.end(), query.size(), true, true, locations, locate);
    bool include_middle = rev_query.size() % 2 == 0;
    return occs + count_one_error_case(rev_csa, rev_query.begin(), rev_query.end(), rev_query.size(), include_middle, false, locations, locate);
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        cout << "Usage " << argv[0] << " text_file [error_num] [max_locations] [post_context] [pre_context]" << endl;
        cout << "    This program constructs a very compact FM-index" << endl;
        cout << "    which supports count, locate, and extract queries." << endl;
        cout << "    text_file      Original text file." << endl;
        cout << "    error_num      The number of errors (mismatches)." << endl;
        cout << "    max_locations  Maximal number of location to report." << endl;
        cout << "    post_context   Maximal length of the reported post-context." << endl;
        cout << "    pre_context    Maximal length of the pre-context." << endl;
        return 1;
    }

    size_t error_num = 0;
    size_t max_locations = 5;
    size_t post_context = 10;
    size_t pre_context = 10;
    if (argc >= 3)
    {
        error_num = atoi(argv[2]);
        if (error_num > 2) // It's already non-negative
        {
            cout << "    error_num must be one of {0, 1, 2}." << endl;
            return 1;
        }
    }
    if (argc >= 4)
    {
        max_locations = atoi(argv[3]);
    }
    if (argc >= 5)
    {
        post_context = atoi(argv[4]);
    }
    if (argc >= 6)
    {
        pre_context = atoi(argv[5]);
    }

    string index_suffix = ".fm9";
    string index_file = string(argv[1]) + index_suffix;
    csa_wt<wt_huff<rrr_vector<127>>, 512, 1024> fm_index;

    if (!load_from_file(fm_index, index_file))
    {
        ifstream in(argv[1]);
        if (!in)
        {
            cout << "ERROR: File " << argv[1] << " does not exist. Exit." << endl;
            return 1;
        }
        cout << "No index " << index_file << " located. Building index now." << endl;
        construct(fm_index, argv[1], 1);     // generate index
        store_to_file(fm_index, index_file); // save it
    }
    cout << "Index construction complete, index requires " << size_in_mega_bytes(fm_index) << " MiB." << endl;

    ifstream index_file_stream(argv[1]);
    string rev_index_string((istreambuf_iterator<char>(index_file_stream)), istreambuf_iterator<char>());
    reverse(rev_index_string.begin(), rev_index_string.end());

    string rev_index_file = "reversed_" + index_file;
    csa_wt<wt_huff<rrr_vector<127>>, 512, 1024> rev_fm_index;

    if (!load_from_file(rev_fm_index, rev_index_file))
    {
        cout << "No reversed index " << index_file << " located. Building reversed index now." << endl;
        construct_im(rev_fm_index, rev_index_string, 1); // generate index
        store_to_file(rev_fm_index, rev_index_file);     // save it
    }

    string temp_query = "bcd";
    // size_t s, l_fwd_res, r_fwd_res, l_bwd_res, r_bwd_res;
    // string::iterator begin = temp.begin();
    // string::iterator end = temp.end();

    // string::iterator begin = temp.begin();
    // size_t count = 0, size = 0;
    // size_t l_fwd_res, r_fwd_res, l_bwd_res, r_bwd_res;
    // bidirectional_search(fm_index, 0, fm_index.size() - 1, 0, fm_index.size() - 1, 'b', l_fwd_res, r_fwd_res, l_bwd_res, r_bwd_res);

    bool do_locate = false;
    int_vector<64> locations(0);
    if (max_locations > 0)
        do_locate = true;

    cout << "Input search terms and press Ctrl-D to exit." << endl;
    string prompt = "\e[0;32m>\e[0m ";
    cout << prompt;
    string query, rev_query;
    while (getline(cin, query))
    {
        size_t m = query.size();

        string rev_query = query;
        reverse(rev_query.begin(), rev_query.end());

        size_t occs;
        switch (error_num)
        {
        case 1:
            occs = handle_one_error(fm_index, rev_fm_index, query, rev_query, locations, do_locate);
            break;
        case 2: // TODO
        default:
            occs = sdsl::count(fm_index, query.begin(), query.end());
            locations = locate(fm_index, query.begin(), query.begin() + m);
        }

        cout << "# of occurrences: " << occs << endl;
        if (occs > 0)
        {
            cout << "Location and context of first occurrences: " << endl;
            // TODO: Uncomment
            // auto locations = locate(fm_index, query.begin(), query.begin() + m);
            sort(locations.begin(), locations.end());
            for (size_t i = 0, pre_extract = pre_context, post_extract = post_context; i < min(occs, max_locations); ++i)
            {
                cout << setw(8) << locations[i] << ": ";
                if (pre_extract > locations[i])
                {
                    pre_extract = locations[i];
                }
                if (locations[i] + m + post_extract > fm_index.size())
                {
                    post_extract = fm_index.size() - locations[i] - m;
                }
                auto s = extract(fm_index, locations[i] - pre_extract, locations[i] + m + post_extract - 1);
                string pre = s.substr(0, pre_extract);
                s = s.substr(pre_extract);
                if (pre.find_last_of('\n') != string::npos)
                {
                    pre = pre.substr(pre.find_last_of('\n') + 1);
                }
                cout << pre;
                cout << "\e[1;31m";
                cout << s.substr(0, m);
                cout << "\e[0m";
                string context = s.substr(m);
                cout << context.substr(0, context.find_first_of('\n')) << endl;
            }
        }
        cout << prompt;
    }

    index_file_stream.close(); // TODO: Check it's fine
    cout << endl;
}
