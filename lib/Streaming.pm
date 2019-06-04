package Streaming;
use Vcf;
use List::MoreUtils qw/all/;


sub create_streamers {
    my ($files,$range,$job_id,$bcftools) = @_;
    
    my @fhs;
    
    foreach ( sort {$a cmp $b } keys %$files ) {
        my %info;
        #fetch only the range from the reference that was assigned to this filehandle.
        my $cmd = "$bcftools query -r $range -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t%INFO/DP\n' " . $files->{$_} . " 2>/dev/null";
        
        #created filehandle from $cmd above directly so we do not have to write to disk
	open (my $OUTPUT,'-|',"$cmd");
        $info{'vcf'} = Vcf->new(fh => $OUTPUT ,_version_set=> '4.0');
        $info{'file'} = $files->{$_};
        push @fhs,\%info;
    }
    
    return sub {
        my ($chrom,$pos) = @_;
        my %line;
        my @result;

        
        foreach my $info( @fhs) {
            while (1) {

                #if we have last entry as EOF, simply return with ''
                if ( $info->{'last'} eq 'EOF') {
                    push @result,{'status'=>'EOF'};
                    last;
                }

                #check to see if we already readed a line but did not process it
                my $data = $info->{'last'} ? $info->{'last'} : $info->{'vcf'}->next_data_array;
                
                $info->{'last'} = undef if $info->{'last'};
                
                #if undef, we are at the end of the 'file'
                if (! $data){
                    $info->{'last'}='EOF';
                    push @result,{'status'=>'EOF'};
                    last;
                }
                else{
                    my $chrom_cur = $data->[0];
                    my $position_cur = $data->[1];
                    
                    # found a position to keep
                    if ($chrom eq $chrom_cur && $pos == $position_cur){

                        my $ref = $data->[3];
                        my $alt = $data->[4];
                        my @status = split /;/ , $data->[6]; #grabbing status
                        
                        
                        push @result,{'ref' => $ref, 'alt' => $alt, 'status' => $status[0]};
                        last;
                    }
                    elsif ( $chrom eq $chrom_cur && $pos < $position_cur) {
                        push @result,{'status' =>''};
                        $info->{'last'}= $data;
                        last;
                    }
                    elsif ( $chrom ne $chrom_cur) {
                        push @result,{'status' =>''};
                        $info->{'last'}= $data;
                        last;
                    }
                }
            }
                
        
            
        }
        
        return @result

    };
    
}



1;
