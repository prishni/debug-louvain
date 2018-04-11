function faid = get_diggs_facet_id(R_type,facet_name)
%    all_dims = [nU,nS,nC,nW,nTP];
    faid = -1;
    switch(R_type)
        case {'R2'} % user-contact (contact)
            if strcmp(facet_name,'user'), faid=1; end
        case {'R23','R234'} % contact,user-story (submit) 
                            % contact,submit,user-story (digg)
            if strcmp(facet_name,'user'), faid=1; end
            if strcmp(facet_name,'story'), faid=2; end
        case {'R123','R1234','R13','R14','R134','R4789','R478','R124'} % story-keyword-topic (topic),contact,submit
                              % topic, contact,submit,digg
            if strcmp(facet_name,'user'), faid=1; end
            if strcmp(facet_name,'story'), faid=2; end
            if strcmp(facet_name,'keyword'), faid=3; end
            if strcmp(facet_name,'topic'), faid=4; end
        case {'R2345','R23456'} % contact,submit,digg,user-story-comment (comment)
                                % contact,submit,digg,comment,user-comment-reply (reply)
            if strcmp(facet_name,'user'), faid=1; end
            if strcmp(facet_name,'story'), faid=2; end
            if strcmp(facet_name,'comment'), faid=3; end
        case {'R12345','R123456','R15','R135','R145','R1345','R125','R5789','R1235','R1245'} % topic,contact,submit,digg,user-story-comment (comment)
                                % topic,contact,submit,digg,comment,user-comment-reply (reply)
            if strcmp(facet_name,'user'), faid=1; end
            if strcmp(facet_name,'story'), faid=2; end
            if strcmp(facet_name,'comment'), faid=3; end
            if strcmp(facet_name,'keyword'), faid=4; end
            if strcmp(facet_name,'topic'), faid=5; end
    end