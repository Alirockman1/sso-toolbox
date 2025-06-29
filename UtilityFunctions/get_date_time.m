function currentTimeStamp = get_date_time
%GET_DATE_TIME Get a time stamp for the current time (yyyy-mm-dd_HH-MM-SS)
%   GET_DATE_TIME checks the clock in the computer to get the current time.
%   The format used for the output is a string in the form: 
%   'Year-Month-Day_Hour-Minute-Second'.
%
%   CURRENTTIMESTAMP = GET_DATE_TIME assigns to CURRENTTIMESTAMP a string
%   with the time-stamp as explained above.
%
%   Output:
%       - CURRENTTIMESTAMP : string 
%
%   See also save_diary_files.
%
%   Copyright 2024 Eduardo Rodrigues Della Noce
%   SPDX-License-Identifier: Apache-2.0

%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
% 
%       http://www.apache.org/licenses/LICENSE-2.0
% 
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.

    c = clock;
    % Year-Month-Day_Hour-Minute-Second
    currentTimeStamp = sprintf('%04g-%02g-%02g_%02g-%02g-%02g',c(1),c(2),c(3),c(4),c(5),round(c(6)));
end