Uses [angel's fork]((https://github.com/angelxuanchang/simple-amt)) of [simple-amt](https://github.com/jcjohnson/simple-amt)
Really the same as the original (just no preview if not qualified and slightly better `get_all_hits.py`) 

`symmetry_object.html` - Template for symmetry verification for objects 

Copy to `simple-amt/examples/symmetry` and run following commands.


## Local testing ##

Generate rendered template

`python render_template.py --html_template=examples/symmetry/symmetry_object.html  --rendered_html=test.html`

Test by running server `python -m SimpleHTTPServer 8085` and going to <localhost:8085/test.html>


## Sandbox ##
1. Launch hits 
    ```
    python launch_hits.py \
      --html_template=examples/symmetry/symmetry_object.html \
      --hit_properties_file=examples/symmetry/symmetry_object.json \
      --input_json_file=examples/symmetry/example_input.jsonl \
      --hit_ids_file=examples/symmetry/hit_ids.txt
    ```
   Test on <https://workersandbox.mturk.com>

1. Check hit progress (shows as 1/1: Number that has 1/1 assignments done, 0/1: Number that is not done)
    
    `python show_hit_progress.py --hit_ids_file=examples/symmetry/hit_ids.txt`  

1. When done, run `get_results.py` to save results
    ```
    python get_results.py \
      --hit_ids_file=examples/symmetry/hit_ids.txt \
      > examples/symmetry/results.txt
    ```

1. Approve all hits

    `python approve_hits.py --hit_ids_file=examples/symmetry/hit_ids.txt`

1. Disable (delete hits) (make sure that that you ran `get_results.py` before this step).  
Seems to also approve outstanding assignments.

    `python disable_hits.py --hit_ids_file=examples/symmetry/hit_ids.txt`


## Production (real turkers) ##
1. Launch hits 
    ```
    python launch_hits.py \
      --prod \
      --html_template=examples/symmetry/symmetry_object.html \
      --hit_properties_file=examples/symmetry/symmetry_object.json \
      --input_json_file=examples/symmetry/input.jsonl \
      --hit_ids_file=examples/symmetry/prod_hit_ids.txt
    ```

1. Check hit progress (shows as 1/1: Number that has 1/1 assignments done, 0/1: Number that is not done)

    `python show_hit_progress.py --prod --hit_ids_file=examples/symmetry/prod_hit_ids.txt`

1. When done, run `get_results.py` to save results
    ```
    python get_results.py \
      --prod \
      --hit_ids_file=examples/symmetry/prod_hit_ids.txt \
      > examples/symmetry/prod_results.txt
    ```

1. Approve all hits
    `python approve_hits.py --prod --hit_ids_file=examples/symmetry/prod_hit_ids.txt`


1. Disable (delete hits) (make sure that that you ran `get_results.py` before this step).  
Seems to also approve outstanding assignments.
    `python disable_hits.py --prod --hit_ids_file=examples/symmetry/prod_hit_ids.txt`


