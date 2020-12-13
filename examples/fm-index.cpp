#include <sdsl/suffix_arrays.hpp>
#include <string>
#include <iostream>
#include <algorithm>
#include <vector>

using namespace sdsl;
using namespace std;

//-----------  1-error -----------//
typedef csa_wt<wt_blcd<>, 32, 32, sa_order_sa_sampling<>, isa_sampling<>, byte_alphabet> t_csa;
template <class t_csa, class t_rac, class t_pat_iter>
typename t_csa::size_type count_one_error_case(const t_csa &csa, typename t_csa::size_type left_window, typename t_csa::size_type right_window,
                                               t_pat_iter begin, t_pat_iter end, bool include_middle, bool case_a, t_rac &locations, bool locate)
{
    typename t_csa::size_type m = end - begin;
    typename t_csa::size_type x = (m + 1) / 2;
    if (!include_middle)
        x--;

    if (end - begin > (typename std::iterator_traits<t_pat_iter>::difference_type)csa.size())
        return 0;

    typename t_csa::char_type curr_char;
    typename t_csa::size_type left_res = 0, right_res = 0, left_err_res = 0, right_err_res = 0, result = 0, i = 0, occs = 0;
    size_t locations_size;
    //find SA of right half P[x...m]
    backward_search(csa, left_window, right_window, begin + x, end, left_res, right_res);

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
                                locations[locations_size + k] = csa.size() - 1 - csa[left_err_res + k] - m;
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

template <class t_csa, class t_rac, class t_pat_iter>
typename t_csa::size_type
handle_one_error(
    const t_csa &csa,
    const t_csa &rev_csa,
    typename t_csa::size_type left_window,
    typename t_csa::size_type right_window,
    t_pat_iter begin,
    t_pat_iter end,
    t_pat_iter rev_begin,
    t_pat_iter rev_end,
    t_rac &locations,
    bool locate)
{
    size_t occs = 0;
    occs = count_one_error_case(csa, left_window, right_window, begin, end, true, true, locations, locate);
    bool include_middle = (rev_end - rev_begin) % 2 == 0;
    return occs + count_one_error_case(rev_csa, left_window, right_window, rev_begin, rev_end, include_middle, false, locations, locate);
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
    return handle_one_error(csa, rev_csa, 0, csa.size() - 1, query.begin(), query.end(), rev_query.begin(), rev_query.end(), locations, locate);
}


//-----------  2-errors -----------//

template <class t_csa, class t_rac, class t_pat_iter>
typename t_csa::size_type count_two_errors_case_a(const t_csa &csa,
                                                t_pat_iter begin,
                                                t_pat_iter end,
                                                t_rac &locations,
                                                bool locate)
{
    size_t locations_size;
    typename t_csa::size_type m = end - begin;
    typename t_csa::size_type s_1 = m / 3;
    typename t_csa::size_type s_2 = m - s_1;

    if (end - begin > (typename std::iterator_traits<t_pat_iter>::difference_type)csa.size())
        return 0;

    typename t_csa::char_type curr_first_char, curr_second_char;
    typename t_csa::size_type left_res = 0, right_res = 0, j = 0, left_err_res = 0, right_err_res = 0,left_err2_res=0, right_err2_res=0, result = 0, occs = 0;

    backward_search(csa, 0, csa.size() - 1, begin + s_2, end, left_res, right_res);
    
    if (left_res <= right_res)
    {
        for (j = s_2; j > 1; j--)
        {
            curr_first_char = (typename t_csa::char_type) * (begin + j - 1);
            //possible replacments of first error:
            for (size_t i = 1; i < csa.sigma; i++)
            {
                if (csa.char2comp[curr_first_char] != i)
                {
                    backward_search(csa, left_res, right_res, csa.comp2char[i],left_err_res, right_err_res);
                    
                    for (size_t k = j-1; k>0 ; k--){
                        curr_second_char = (typename t_csa::char_type) * (begin + k-1);
                                    
                        //possible replacments of second error:
                        for (size_t l = 1; l < csa.sigma; l++){
                            if (csa.char2comp[curr_second_char] != l){
                                
                                backward_search(csa, left_err_res, right_err_res,  csa.comp2char[l],left_err2_res, right_err2_res); 
                                occs = backward_search(csa, left_err2_res, right_err2_res, begin, begin+k-1, left_err2_res,right_err2_res);
                                result+= occs;
                                
                                if (locate && occs > 0)
                                {
                                    locations_size = locations.size();
                                    locations.resize(locations_size + occs);
                                    for (typename t_csa::size_type t = 0; t < occs; t++)
                                    {
                                        locations[locations_size + t] = csa[left_err2_res + t];
                                    }
                                }
                            }
                        }
                        backward_search(csa, left_err_res,right_err_res, curr_second_char, left_err_res,right_err_res);
                    }
                }
            }
            backward_search(csa, left_res, right_res, curr_first_char, left_res, right_res);
        }
    }
    return result;
}

template <class t_csa, class t_rac, class t_pat_iter>
typename t_csa::size_type
count_two_errors_case_b(const t_csa &rev_csa,
                        t_pat_iter rev_begin,
                        t_pat_iter rev_end,
                        t_rac &locations,
                        bool locate)
{
    typename t_csa::size_type m = rev_end - rev_begin;
    typename t_csa::size_type s_1 = m / 3; // Original, not the reversed
    typename t_csa::size_type s_2 = m - s_1; // Original, not the reversed

    if (rev_end - rev_begin > (typename std::iterator_traits<t_pat_iter>::difference_type)rev_csa.size())
        return 0;

    typename t_csa::char_type curr_char, curr_char2;
    typename t_csa::size_type rev_left_res = 0, rev_right_res = 0,
                              rev_left_err_res = 0, rev_right_err_res = 0,
                              rev_left_err2_res = 0, rev_right_err2_res = 0,
                              occs = 0, result = 0, locations_size, i, j, k, l, n;

    backward_search(rev_csa, 0, rev_csa.size() - 1, rev_begin + m - s_2, rev_end, rev_left_res, rev_right_res);
    
    if (rev_left_res > rev_right_res)
        return 0;

    for (i = m - s_2 - 1; i > 0; i--)
    {
        curr_char = (typename t_csa::char_type) * (rev_begin + i);
        for (j = 1; j < rev_csa.sigma; j++)
        {
            if (rev_csa.char2comp[curr_char] != j)
            {
                occs = backward_search(rev_csa, rev_left_res, rev_right_res, rev_csa.comp2char[j], rev_left_err_res, rev_right_err_res);

                if (occs == 0)
                    continue;

                for (k = i; k > 0; k--)
                {
                    curr_char2 = (typename t_csa::char_type) * (rev_begin + k - 1);
                    for (l = 1; l < rev_csa.sigma; l++)
                    {
                        if (rev_csa.char2comp[curr_char2] != l)
                        {
                            backward_search(rev_csa, rev_left_err_res, rev_right_err_res, rev_csa.comp2char[l], rev_left_err2_res, rev_right_err2_res);

                            occs = backward_search(rev_csa, rev_left_err2_res, rev_right_err2_res, rev_begin, rev_begin + k - 1, rev_left_err2_res, rev_right_err2_res);
                            
                            result += occs;

                            if (locate && occs > 0)
                            {
                                locations_size = locations.size();
                                locations.resize(locations_size + occs);
                                for (n = 0; n < occs; n++)
                                {
                                    locations[locations_size + n] = rev_csa.size() - 1 - rev_csa[rev_left_err2_res + n] - m;
                                }
                            }
                        }
                    }
                    backward_search(rev_csa, rev_left_err_res, rev_right_err_res, curr_char2, rev_left_err_res, rev_right_err_res);
                }
            }
        }
        backward_search(rev_csa, rev_left_res, rev_right_res, curr_char, rev_left_res, rev_right_res);
    }
    return result;
}

template <class t_csa, class t_rac, class t_pat_iter>
typename t_csa::size_type
count_two_errors_case_c(const t_csa &rev_csa,
                        t_pat_iter rev_begin,
                        t_pat_iter rev_end,
                        t_rac &locations,
                        bool locate)
{
    typename t_csa::size_type m = rev_end - rev_begin;
    typename t_csa::size_type s_1 = m / 3; // Original, not the reversed
    typename t_csa::size_type s_2 = m - s_1; // Original, not the reversed

    if (rev_end - rev_begin > (typename std::iterator_traits<t_pat_iter>::difference_type)rev_csa.size())
        return 0;

    typename t_csa::char_type curr_char, curr_char2;
    typename t_csa::size_type rev_left_res = 0, rev_right_res = 0,
                              rev_left_err_res = 0, rev_right_err_res = 0,
                              rev_left_err2_res = 0, rev_right_err2_res = 0,
                              occs = 0, result = 0, locations_size, i, j, k, l, n;

    backward_search(rev_csa, 0, rev_csa.size() - 1, rev_begin + m - s_1, rev_end, rev_left_res, rev_right_res);
    
    if (rev_left_res > rev_right_res)
        return 0;

    for (i = m - 1 - s_1; i > m - 1 - s_2; i--)
    {
        curr_char = (typename t_csa::char_type) * (rev_begin + i);
        for (j = 1; j < rev_csa.sigma; j++)
        {
            if (rev_csa.char2comp[curr_char] != j)
            {
                backward_search(rev_csa, rev_left_res, rev_right_res, rev_csa.comp2char[j], rev_left_err_res, rev_right_err_res);
                
                occs = backward_search(rev_csa, rev_left_err_res, rev_right_err_res, rev_begin + m - s_2, rev_begin + i, rev_left_err_res, rev_right_err_res);
                
                if (occs == 0)
                    continue;

                for (k = m - s_2; k > 0; k--)
                {
                    curr_char2 = (typename t_csa::char_type) * (rev_begin + k - 1);
                    for (l = 1; l < rev_csa.sigma; l++)
                    {
                        if (rev_csa.char2comp[curr_char2] != l)
                        {
                            backward_search(rev_csa, rev_left_err_res, rev_right_err_res, rev_csa.comp2char[l], rev_left_err2_res, rev_right_err2_res);

                            occs = backward_search(rev_csa, rev_left_err2_res, rev_right_err2_res, rev_begin, rev_begin + k - 1, rev_left_err2_res, rev_right_err2_res);
                            
                            result += occs;

                            if (locate && occs > 0)
                            {
                                locations_size = locations.size();
                                locations.resize(locations_size + occs);
                                for (n = 0; n < occs; n++)
                                {
                                    locations[locations_size + n] = rev_csa.size() - 1 - rev_csa[rev_left_err2_res + n] - m;
                                }
                            }
                        }
                    }
                    backward_search(rev_csa, rev_left_err_res, rev_right_err_res, curr_char2, rev_left_err_res, rev_right_err_res);
                }
            }
        }
        backward_search(rev_csa, rev_left_res, rev_right_res, curr_char, rev_left_res, rev_right_res);
    }
    return result;
}

template <class t_csa, class t_rac, class t_pat_iter>
typename t_csa::size_type
count_two_errors_case_d(const t_csa &csa,
                        const t_csa &rev_csa,
                        t_pat_iter begin,
                        t_pat_iter end,
                        t_rac &locations,
                        bool locate)
{
    typename t_csa::size_type m = end - begin;
    typename t_csa::size_type s_1 = m / 3;
    typename t_csa::size_type s_2 = m - s_1;

    if (end - begin > (typename std::iterator_traits<t_pat_iter>::difference_type)csa.size())
        return 0;

    typename t_csa::char_type curr_char, curr_char2;
    typename t_csa::size_type left_res = 0, right_res = 0, rev_left_res = 0, rev_right_res = 0,
                              left_err_res = 0, right_err_res = 0, rev_left_err_res = 0, rev_right_err_res = 0,
                              left_err2_res = 0, right_err2_res = 0, rev_left_err2_res = 0, rev_right_err2_res = 0,
                              occs = 0, result = 0, locations_size, i, j, k, l, n;

    // First obtain the the SA range of P[s1+1..s2] and then the SA’ range using forward search.
    bidirectional_search_forward(csa, rev_csa, 0, csa.size() - 1, 0, rev_csa.size() - 1, begin + s_1, begin + s_2,
                                 left_res, right_res, rev_left_res, rev_right_res);

    if (left_res > right_res)
        return 0;

    // For each i=s1-1,...,0, we apply backward search to compute the SA range of P[1..i−1]e1P[i+1..s2-1].
    for (i = s_1; i > 0; i--)
    {
        curr_char = (typename t_csa::char_type) * (begin + i - 1);
        for (j = 1; j < csa.sigma; j++)
        {
            if (csa.char2comp[curr_char] != j)
            {
                bidirectional_search(csa, left_res, right_res, rev_left_res, rev_right_res, csa.comp2char[j],
                                     left_err_res, right_err_res, rev_left_err_res, rev_right_err_res);

                occs = bidirectional_search_backward(csa, rev_csa, left_err_res, right_err_res, rev_left_err_res, rev_right_err_res, begin, begin + i - 1,
                                                     left_err_res, right_err_res, rev_left_err_res, rev_right_err_res);

                if (occs == 0)
                    continue;

                // For each k=s2,...,m, we apply forward search to compute the SA range of P[1..i−1]e1P[i..j−1]e2P[j..m] for all possible e2.
                for (k = s_2; k < m; k++)
                {
                    curr_char2 = (typename t_csa::char_type) * (begin + k);
                    for (l = 1; l < csa.sigma; l++)
                    {
                        if (csa.char2comp[curr_char2] != l)
                        {
                            bidirectional_search(rev_csa, rev_left_err_res, rev_right_err_res, left_err_res, right_err_res, csa.comp2char[l],
                                                 rev_left_err2_res, rev_right_err2_res, left_err2_res, right_err2_res);

                            occs = bidirectional_search_forward(csa, rev_csa, left_err2_res, right_err2_res, rev_left_err2_res, rev_right_err2_res, begin + k + 1, end,
                                                                left_err2_res, right_err2_res, rev_left_err2_res, rev_right_err2_res);

                            result += occs;

                            if (locate && occs > 0)
                            {
                                locations_size = locations.size();
                                locations.resize(locations_size + occs);
                                for (n = 0; n < occs; n++)
                                    locations[locations_size + n] = csa[left_err2_res + n];
                            }
                        }
                    }
                    bidirectional_search(rev_csa, rev_left_err_res, rev_right_err_res, left_err_res, right_err_res, curr_char2,
                                         rev_left_err_res, rev_right_err_res, left_err_res, right_err_res);
                }
            }
        }
        //return char at i index in P to the original one
        bidirectional_search(csa, left_res, right_res, rev_left_res, rev_right_res, curr_char,
                             left_res, right_res, rev_left_res, rev_right_res);
    }
    return result;
}

template <class t_csa, class t_rac>
typename t_csa::size_type
handle_two_errors(
    const t_csa &csa,
    const t_csa &rev_csa,
    string query,
    string rev_query,
    t_rac &locations,
    bool locate)
{
    size_t occs = 0;
    occs += count_two_errors_case_a(csa, query.begin(), query.end(), locations, locate);
    occs += count_two_errors_case_b(rev_csa, rev_query.begin(), rev_query.end(), locations, locate);
    occs += count_two_errors_case_c(rev_csa, rev_query.begin(), rev_query.end(), locations, locate);
    occs += count_two_errors_case_d(csa, rev_csa, query.begin(), query.end(), locations, locate);
    return occs;
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
    t_csa fm_index;

    if (!load_from_file(fm_index, index_file))
    {
        ifstream in(argv[1]);
        if (!in)
        {
            cout << "ERROR: File " << argv[1] << " does not exist. Exit." << endl;
            return 1;
        }
        construct(fm_index, argv[1], 1);     // generate index
        store_to_file(fm_index, index_file); // save it
    }

    ifstream index_file_stream(argv[1]);
    string rev_index_string((istreambuf_iterator<char>(index_file_stream)), istreambuf_iterator<char>());
    reverse(rev_index_string.begin(), rev_index_string.end());

    string rev_index_file = "reversed_" + index_file;
    t_csa rev_fm_index;

    if (!load_from_file(rev_fm_index, rev_index_file))
    {
        construct_im(rev_fm_index, rev_index_string, 1); // generate index
        store_to_file(rev_fm_index, rev_index_file);     // save it
    }

    bool do_locate = false;
    if (max_locations > 0)
        do_locate = true;

    string query, rev_query;
    while (getline(cin, query))
    {
        size_t m = query.size();
        int_vector<64> locations(0);
        string rev_query = query;
        reverse(rev_query.begin(), rev_query.end());

        size_t occs;
        switch (error_num)
        {
        case 1:
            occs = handle_one_error(fm_index, rev_fm_index, query, rev_query, locations, do_locate);
            break;
        case 2:
            occs = handle_two_errors(fm_index, rev_fm_index, query, rev_query, locations, do_locate);
            break;
        default:
            occs = sdsl::count(fm_index, query.begin(), query.end());
            locations = locate(fm_index, query.begin(), query.begin() + m);
        }

        cout << query << " : " << occs << endl;
        if (occs > 0)
        {
            if (do_locate)
                cout << "Location and context of first occurrences:" << endl;
            sort(locations.begin(), locations.end());
            for (size_t i = 0, pre_extract = pre_context, post_extract = post_context; i < min(occs, max_locations); ++i)
            {
                cout << "  " << locations[i] << ": ";
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
                cout << s.substr(0, m);
                cout << endl;
            }
        }
    }

    index_file_stream.close(); 
}
