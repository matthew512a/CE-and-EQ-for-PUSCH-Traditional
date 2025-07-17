function channel = getChannelInfo(profile,fd,tm,Fs)

% Get channel config
DelayProfile                        = profile;
MaximumDopplerShift                 = fd;
DelaySpread                         = tm;
SampleRate                          = Fs;
NumTransmitAntennas                 = 1;
NumReceiveAntennas                  = 1;
InitailTime                         = 0;

% Get channel parameters
channel                             = nrTDLChannel;
channel.DelayProfile                = DelayProfile;
channel.DelaySpread                 = DelaySpread;
channel.MaximumDopplerShift         = MaximumDopplerShift;
channel.SampleRate                  = SampleRate;
channel.MIMOCorrelation             = 'Low';
channel.Polarization                = "Cross-Polar";
channel.TransmissionDirection       = "Uplink";
channel.NumTransmitAntennas         = NumTransmitAntennas;
channel.NumReceiveAntennas          = NumReceiveAntennas;
channel.NormalizePathGains          = true;
channel.InitialTime                 = InitailTime;
channel.NumSinusoids                = 48;
channel.RandomStream                = 'mt19937ar with seed';
channel.Seed                        = 73;
channel.NormalizeChannelOutputs     = true;
end