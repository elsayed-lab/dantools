use Test::More;
eval "use Test::Pod 1.00";
plan skip_all => "Test::Pod 1.00 required for testing the validity of the POD." if $@;
all_pod_files_ok();
